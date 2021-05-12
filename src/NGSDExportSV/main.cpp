#include "ToolBase.h"
#include "NGSD.h"
#include "Exceptions.h"
#include "Helper.h"
#include "BedpeFile.h"

class ConcreteTool
		: public ToolBase
{
	Q_OBJECT

public:
	ConcreteTool(int& argc, char *argv[])
		: ToolBase(argc, argv)
	{
	}

	virtual void setup()
	{
		setDescription("Exports structural variants of the NGSD into BEDPE file(s).");
		addOutfile("out", "Output BEDPE file prefix for structural variants.", false);

		//optional
//		addInt("minCount", "Minimal NGSD count for a SV to be exported.", true, 0);
//		addFloat("minAF", "Minimal NGSD allele frequency for a SV to be exported.", true, 0.0);
		addString("system", "Processing system which should be exported.", true, "TruSeqPCRfree");
		addFlag("test", "Uses the test database instead of on the production database.");
		addFlag("debug", "Provide additional information in STDOUT (e.g. query runtime)");
		addString("types", "Comma seperated list of SV types which should be exported", true, "DEL,DUP,INV,INS,BND");

		changeLog(2021, 5, 11, "Initial version.");
	}

	virtual void main()
	{
		//init
		NGSD db(getFlag("test"));
		QTextStream out(stdout);
		QTime timer;
		timer.start();
		bool debug = getFlag("debug");
		QTime init_timer, del_timer, dup_timer, inv_timer, ins_timer, bnd_timer;
		int time_init=0, time_sum_del=0, time_sum_dup=0, time_sum_inv=0, time_sum_ins=0, time_sum_bnd=0;
		int n_del=0, n_dup=0, n_inv=0, n_ins=0, n_bnd=0;
		int n_del_exported=0, n_dup_exported=0, n_inv_exported=0, n_ins_exported=0, n_bnd_exported=0;
		int sample_count;

		//parse types to export
		QStringList type_strings = getString("types").split(',');
		QList<StructuralVariantType> types;
		foreach (const QString& type_string, type_strings)
		{
			types << StructuralVariantTypeFromString(type_string.trimmed());
		}


		BedpeFile bedpe_file_dummy; //dummy file to recieve svs
		bedpe_file_dummy.setAnnotationHeaders(QList<QByteArray>() << "TYPE" << "SYSTEM" << "NGSD_COUNT" << "NGSD_AF");

		if (debug) init_timer.start();

		// create temporary tables for each SV type filtered by the current processing system
		QByteArray table_prefix;

		// create temp dbs in memory
		QByteArray db_engine; // = "ENGINE=MEMORY ";

		// get processing system of current sample
		int processing_system_id = db.processingSystemId(getString("system"));
		QByteArray processing_system = db.getProcessingSystemData(processing_system_id).name_short.toUtf8();

		// create a temp table with all valid ps ids ignoring bad/merged samples
		SqlQuery temp_table = db.getQuery();
		temp_table.exec("CREATE TEMPORARY TABLE temp_valid_sv_cs_ids " + db_engine + "SELECT sc.id FROM sv_callset sc INNER JOIN processed_sample ps ON sc.processed_sample_id = ps.id WHERE ps.processing_system_id = " + QByteArray::number(processing_system_id) + " AND ps.quality != 'bad' AND NOT EXISTS (SELECT 1 FROM merged_processed_samples mps WHERE mps.processed_sample_id = sc.processed_sample_id)");

		// get number of valid callset ids (= known samples)
		sample_count = db.getValue("SELECT COUNT(*) FROM temp_valid_sv_cs_ids").toInt();

		// generate joined tables of all SV types which were called on the same processing system.d
		SqlQuery create_temp_table = db.getQuery();
		QByteArrayList table_names;
		table_names << "sv_deletion" << "sv_duplication" << "sv_inversion";
		foreach (QByteArray table_name, table_names)
		{
			// DEL, DUP, INV
			create_temp_table.exec("CREATE TEMPORARY TABLE temp_" + table_name + " " + db_engine + "SELECT sv.id, sv.sv_callset_id, sv.chr, sv.start_min, sv.start_max, sv.end_min, sv.end_max FROM "
								   + table_name + " sv INNER JOIN temp_valid_sv_cs_ids tt ON sv.sv_callset_id = tt.id");
		}
		//INS
		create_temp_table.exec("CREATE TEMPORARY TABLE temp_sv_insertion " + db_engine + "SELECT sv.id, sv.sv_callset_id, sv.chr, sv.pos, sv.ci_upper FROM sv_insertion sv INNER JOIN temp_valid_sv_cs_ids tt ON sv.sv_callset_id = tt.id");
		//BND
		create_temp_table.exec("CREATE TEMPORARY TABLE temp_sv_translocation " + db_engine + "SELECT sv.id, sv.sv_callset_id, sv.chr1, sv.start1, sv.end1, sv.chr2, sv.start2, sv.end2 FROM sv_translocation sv INNER JOIN temp_valid_sv_cs_ids tt ON sv.sv_callset_id = tt.id");

		// create indices for exact and overlap matching
		SqlQuery create_index = db.getQuery();
		foreach (QByteArray table_name, table_names)
		{
			// DEL, DUP, INV
			create_index.exec("CREATE INDEX `exact_match` ON temp_" + table_name + "(`chr`, `start_min`, `start_max`, `end_min`, `end_max`)");
		}
		//INS
		create_index.exec("CREATE INDEX `match` ON temp_sv_insertion(`chr`, `pos`, `ci_upper`)");
		//BND
		create_index.exec("CREATE INDEX `match` ON temp_sv_translocation(`chr1`, `start1`, `end1`, `chr2`, `start2`, `end2`)");

		// set prefix for temp tables
		table_prefix = "temp_";



		//prepare queries
		QByteArray select_count = "SELECT id FROM ";


		QByteArray exact_match_del_dup_inv = "WHERE sv.chr = :0 AND sv.start_min <= :1 AND :2 <= sv.start_max AND sv.end_min <= :3 AND :4 <= sv.end_max";

		SqlQuery count_sv_deletion_em = db.getQuery();
		count_sv_deletion_em.prepare(select_count + table_prefix + "sv_deletion sv " + exact_match_del_dup_inv);
		SqlQuery count_sv_duplication_em = db.getQuery();
		count_sv_duplication_em.prepare(select_count + table_prefix + "sv_duplication sv " + exact_match_del_dup_inv);
		SqlQuery count_sv_inversion_em = db.getQuery();
		count_sv_inversion_em.prepare(select_count + table_prefix + "sv_inversion sv " + exact_match_del_dup_inv);

		QByteArray match_ins = table_prefix + "sv_insertion sv WHERE sv.chr = :0 AND sv.pos <= :1 AND :2 <= (sv.pos + sv.ci_upper)";

		SqlQuery count_sv_insertion_m = db.getQuery();
		count_sv_insertion_m.prepare(select_count + match_ins);

		QByteArray match_bnd = table_prefix + "sv_translocation sv WHERE sv.chr1 = :0 AND sv.start1 <= :1 AND :2 <= sv.end1 AND sv.chr2 = :3 AND sv.start2 <= :4 AND :5 <= sv.end2";

		SqlQuery count_sv_translocation_em = db.getQuery();
		count_sv_translocation_em.prepare(select_count + match_bnd);

		if (debug)
		{
			time_init = init_timer.elapsed();
			out << "\t init: " << (double) time_init/1000.00 << "s " << endl;
		}

		if (types.contains(StructuralVariantType::DEL))
		{
			/****** deletions ********/
			if (debug) del_timer.start();
			// open output file
			QSharedPointer<QFile> output_file = Helper::openFileForWriting(getOutfile("out") + "_DEL.bedpe",false,false);
			output_file->write("#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + bedpe_file_dummy.annotationHeaders().join("\t") + "\n");

			QStringList sv_ids = db.getValues("SELECT id FROM " + table_prefix + "sv_deletion");
			if (debug) out << sv_ids.count() << " deletions to parse" << endl;
			QSet<int> processed_ids;
			int lines_in_buffer = 0;

			foreach (const QString& id_string, sv_ids)
			{
				int id = Helper::toInt(id_string);
				int ngsd_count_em = 0;

				// count SVs
				n_del++;

				if (debug && (n_del%1000 == 0)) out << n_del << " deletion parsed" << endl;

				// skip processed ids
				if (processed_ids.contains(id))
				{
					//remove id to keep set small
					processed_ids.remove(id);
					continue;
				}

				//get SV
				BedpeLine sv = db.structuralVariant(id, StructuralVariantType::DEL, bedpe_file_dummy, true);

				// get exact matches
				SqlQuery query_em = count_sv_deletion_em;

				// bind values
				query_em.bindValue(0, sv.chr1().strNormalized(true));
				query_em.bindValue(1, sv.end1());
				query_em.bindValue(2, sv.start1());
				query_em.bindValue(3, sv.end2());
				query_em.bindValue(4, sv.start2());

				// execute query
				query_em.exec();

				// parse result
				while(query_em.next())
				{
					ngsd_count_em++;
					processed_ids << query_em.value(0).toInt();
				}

				// write to file
				QList<QByteArray> annotations;
				annotations << "DEL";
				annotations << processing_system;
				annotations << QByteArray::number(ngsd_count_em);
				if (sample_count > 0)
				{
					annotations << QByteArray::number((double) ngsd_count_em / (double) sample_count);
				}
				else
				{
					annotations << "0.000";
				}
				sv.setAnnotations(QList<QByteArray>());
				output_file->write(sv.toTsv() + "\t" + annotations.join('\t') + "\n");
				n_del_exported++;
				lines_in_buffer++;

				if(lines_in_buffer%1000 == 0)
				{
					output_file->flush();
					lines_in_buffer = 0;
				}

			}

			if (debug) time_sum_del += del_timer.elapsed();

			output_file->close();
		}

		if (types.contains(StructuralVariantType::DUP))
		{
			/****** duplication ********/
			if (debug) dup_timer.start();
			// open output file
			QSharedPointer<QFile> output_file = Helper::openFileForWriting(getOutfile("out") + "_DUP.bedpe",false,false);
			output_file->write("#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + bedpe_file_dummy.annotationHeaders().join("\t") + "\n");

			QStringList sv_ids = db.getValues("SELECT id FROM " + table_prefix + "sv_duplication");
			if (debug) out << sv_ids.count() << " duplications to parse" << endl;
			QSet<int> processed_ids;
			int lines_in_buffer = 0;

			foreach (const QString& id_string, sv_ids)
			{
				int id = Helper::toInt(id_string);
				int ngsd_count_em = 0;

				// count SVs
				n_dup++;

				if (debug && (n_dup%1000 == 0)) out << n_dup << " duplication parsed" << endl;

				// skip processed ids
				if (processed_ids.contains(id))
				{
					//remove id to keep set small
					processed_ids.remove(id);
					continue;
				}

				//get SV
				BedpeLine sv = db.structuralVariant(id, StructuralVariantType::DUP, bedpe_file_dummy, true);

				// get exact matches
				SqlQuery query_em = count_sv_duplication_em;

				// bind values
				query_em.bindValue(0, sv.chr1().strNormalized(true));
				query_em.bindValue(1, sv.end1());
				query_em.bindValue(2, sv.start1());
				query_em.bindValue(3, sv.end2());
				query_em.bindValue(4, sv.start2());

				// execute query
				query_em.exec();

				// parse result
				while(query_em.next())
				{
					ngsd_count_em++;
					processed_ids << query_em.value(0).toInt();
				}

				// write to file
				QList<QByteArray> annotations;
				annotations << "DUP";
				annotations << processing_system;
				annotations << QByteArray::number(ngsd_count_em);
				if (sample_count > 0)
				{
					annotations << QByteArray::number((double) ngsd_count_em / (double) sample_count);
				}
				else
				{
					annotations << "0.000";
				}
				sv.setAnnotations(QList<QByteArray>());
				output_file->write(sv.toTsv() + "\t" + annotations.join('\t') + "\n");
				n_dup_exported++;
				lines_in_buffer++;

				if(lines_in_buffer%1000 == 0)
				{
					output_file->flush();
					lines_in_buffer = 0;
				}

			}

			if (debug) time_sum_dup += dup_timer.elapsed();

			output_file->close();
		}


		if (types.contains(StructuralVariantType::INV))
		{
			/****** inversion ********/
			if (debug) inv_timer.start();
			// open output file
			QSharedPointer<QFile> output_file = Helper::openFileForWriting(getOutfile("out") + "_INV.bedpe",false,false);
			output_file->write("#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + bedpe_file_dummy.annotationHeaders().join("\t") + "\n");

			QStringList sv_ids = db.getValues("SELECT id FROM " + table_prefix + "sv_inversion");
			if (debug) out << sv_ids.count() << " inversions to parse" << endl;
			QSet<int> processed_ids;
			int lines_in_buffer = 0;

			foreach (const QString& id_string, sv_ids)
			{
				int id = Helper::toInt(id_string);
				int ngsd_count_em = 0;

				// count SVs
				n_inv++;

				if (debug && (n_inv%1000 == 0)) out << n_inv << " inversions parsed" << endl;

				// skip processed ids
				if (processed_ids.contains(id))
				{
					//remove id to keep set small
					processed_ids.remove(id);
					continue;
				}

				//get SV
				BedpeLine sv = db.structuralVariant(id, StructuralVariantType::INV, bedpe_file_dummy, true);

				// get exact matches
				SqlQuery query_em = count_sv_inversion_em;

				// bind values
				query_em.bindValue(0, sv.chr1().strNormalized(true));
				query_em.bindValue(1, sv.end1());
				query_em.bindValue(2, sv.start1());
				query_em.bindValue(3, sv.end2());
				query_em.bindValue(4, sv.start2());

				// execute query
				query_em.exec();

				// parse result
				while(query_em.next())
				{
					ngsd_count_em++;
					processed_ids << query_em.value(0).toInt();
				}

				// write to file
				QList<QByteArray> annotations;
				annotations << "INV";
				annotations << processing_system;
				annotations << QByteArray::number(ngsd_count_em);
				if (sample_count > 0)
				{
					annotations << QByteArray::number((double) ngsd_count_em / (double) sample_count);
				}
				else
				{
					annotations << "0.000";
				}
				sv.setAnnotations(QList<QByteArray>());
				output_file->write(sv.toTsv() + "\t" + annotations.join('\t') + "\n");
				n_inv_exported++;
				lines_in_buffer++;

				if(lines_in_buffer%1000 == 0)
				{
					output_file->flush();
					lines_in_buffer = 0;
				}
			}

			if (debug) time_sum_inv += inv_timer.elapsed();

			output_file->close();
		}


		if (types.contains(StructuralVariantType::INS))
		{
			/****** insertion ********/
			if (debug) ins_timer.start();
			// open output file
			QSharedPointer<QFile> output_file = Helper::openFileForWriting(getOutfile("out") + "_INS.bedpe",false,false);
			output_file->write("#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + bedpe_file_dummy.annotationHeaders().join("\t") + "\n");

			QStringList sv_ids = db.getValues("SELECT id FROM " + table_prefix + "sv_insertion");
			if (debug) out << sv_ids.count() << " insertions to parse" << endl;
			QSet<int> processed_ids;
			int lines_in_buffer = 0;

			foreach (const QString& id_string, sv_ids)
			{
				int id = Helper::toInt(id_string);
				int ngsd_count_em = 0;

				// count SVs
				n_ins++;

				if (debug && (n_ins%1000 == 0)) out << n_ins << " insertions parsed" << endl;

				// skip processed ids
				if (processed_ids.contains(id))
				{
					//remove id to keep set small
					processed_ids.remove(id);
					continue;
				}

				//get SV
				BedpeLine sv = db.structuralVariant(id, StructuralVariantType::INS, bedpe_file_dummy, true);

				//get min and max position
				int min_pos = std::min(sv.start1(), sv.start2());
				int max_pos = std::max(sv.end1(), sv.end2());
				//get exact matches
				// bind values
				count_sv_insertion_m.bindValue(0, sv.chr1().strNormalized(true));
				count_sv_insertion_m.bindValue(1, max_pos);
				count_sv_insertion_m.bindValue(2, min_pos);

				// execute query
				count_sv_insertion_m.exec();

				// parse result
				while(count_sv_insertion_m.next())
				{
					ngsd_count_em++;
					processed_ids << count_sv_insertion_m.value(0).toInt();
				}

				// write to file
				QList<QByteArray> annotations;
				annotations << "INS";
				annotations << processing_system;
				annotations << QByteArray::number(ngsd_count_em);
				if (sample_count > 0)
				{
					annotations << QByteArray::number((double) ngsd_count_em / (double) sample_count);
				}
				else
				{
					annotations << "0.000";
				}
				sv.setAnnotations(QList<QByteArray>());
				output_file->write(sv.toTsv() + "\t" + annotations.join('\t') + "\n");
				n_ins_exported++;
				lines_in_buffer++;

				if(lines_in_buffer%1000 == 0)
				{
					output_file->flush();
					lines_in_buffer = 0;
				}
			}

			if (debug) time_sum_ins += ins_timer.elapsed();

			output_file->close();
		}




		if (types.contains(StructuralVariantType::BND))
		{
			/****** translocations ********/
			if (debug) bnd_timer.start();
			// open output file
			QSharedPointer<QFile>  output_file = Helper::openFileForWriting(getOutfile("out") + "_BND.bedpe",false,false);
			output_file->write("#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + bedpe_file_dummy.annotationHeaders().join("\t") + "\n");

			QStringList sv_ids = db.getValues("SELECT id FROM " + table_prefix + "sv_translocation");
			if (debug) out << sv_ids.count() << " translocations to parse" << endl;
			QSet<int> processed_ids;
			int lines_in_buffer = 0;

			foreach (const QString& id_string, sv_ids)
			{
				int id = Helper::toInt(id_string);
				int ngsd_count_em = 0;

				// count SVs
				n_bnd++;

				if (debug && (n_bnd%1000 == 0)) out << n_bnd << " translocations parsed" << endl;

				// skip processed ids
				if (processed_ids.contains(id))
				{
					//remove id to keep set small
					processed_ids.remove(id);
					continue;
				}

				//get SV
				BedpeLine sv = db.structuralVariant(id, StructuralVariantType::BND, bedpe_file_dummy, true);

				count_sv_translocation_em.bindValue(0, sv.chr1().strNormalized(true));
				count_sv_translocation_em.bindValue(1, sv.end1());
				count_sv_translocation_em.bindValue(2, sv.start1());
				count_sv_translocation_em.bindValue(3, sv.chr2().strNormalized(true));
				count_sv_translocation_em.bindValue(4, sv.end2());
				count_sv_translocation_em.bindValue(5, sv.start2());
				// execute query
				count_sv_translocation_em.exec();

				// parse result
				while(count_sv_translocation_em.next())
				{
					ngsd_count_em++;
					processed_ids << count_sv_translocation_em.value(0).toInt();
				}

				// write to file
				QList<QByteArray> annotations;
				annotations << "BND";
				annotations << processing_system;
				annotations << QByteArray::number(ngsd_count_em);
				if (sample_count > 0)
				{
					annotations << QByteArray::number((double) ngsd_count_em / (double) sample_count);
				}
				else
				{
					annotations << "0.000";
				}
				sv.setAnnotations(QList<QByteArray>());
				output_file->write(sv.toTsv() + "\t" + annotations.join('\t') + "\n");
				n_bnd_exported++;
				lines_in_buffer++;

				if(lines_in_buffer%1000 == 0)
				{
					output_file->flush();
					lines_in_buffer = 0;
				}

			}

			if (debug) time_sum_bnd += bnd_timer.elapsed();

			output_file->close();
		}



		// debug summary
		if (debug)
		{
			out << "Debug-Output: \n SQL query runtime / SVs:\n";
			out << "\t init: " << (double) time_init/1000.00 << "s \n";
			out << "\t DEL: " << (double) time_sum_del/1000.00 << "s / " << n_del << " (" << n_del_exported << " exported )\n";
			out << "\t DUP: " << (double) time_sum_dup/1000.00 << "s / " << n_dup << " (" << n_dup_exported << " exported )\n";
			out << "\t INV: " << (double) time_sum_inv/1000.00 << "s / " << n_inv << " (" << n_inv_exported << " exported )\n";
			out << "\t INS: " << (double) time_sum_ins/1000.00 << "s / " << n_ins << " (" << n_ins_exported << " exported )\n";
			out << "\t BND: " << (double) time_sum_bnd/1000.00 << "s / " << n_bnd << " (" << n_bnd_exported << " exported )\n";
		}

		out << "finished (" << Helper::elapsedTime(timer) << ") " << endl;
	}
};

#include "main.moc"

int main(int argc, char *argv[])
{
	ConcreteTool tool(argc, argv);
	return tool.execute();
}
