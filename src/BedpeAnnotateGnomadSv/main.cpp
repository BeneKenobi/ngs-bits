#include "BedFile.h"
#include "ToolBase.h"
#include "NGSHelper.h"
#include "Settings.h"
#include "VcfFile.h"
#include "BedpeFile.h"
#include <QTextStream>
#include <QFileInfo>
#include <QElapsedTimer>

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
		setDescription("Annotates a BEDPE file with information from an Gnomad BED file.");
		addInfile("bed", "Gnomad BED file that is used as annotation source.", false);
		//optional
		addInfile("in", "Input BEDPE file. If unset, reads from STDIN.", true);
		addOutfile("out", "Output BEDPE file. If unset, writes to STDOUT.", true);
		addInt("offset", "Distance for fuzzy match of SVs.", true, 0);

		changeLog(2020, 3, 20, "Initial commit.");
	}

	virtual void main()
	{
		//init
		QString in = getInfile("in");
		QString bed = getInfile("bed");
		QString out = getOutfile("out");
		int offset = getInt("offset");
		if (offset < 0) THROW(ArgumentException, "Offset hast to be positive!");



		//load annoation database
		BedFile gnomad_file;
		gnomad_file.load(bed);
		if (!gnomad_file.isSorted()) gnomad_file.sort();
		ChromosomalIndex<BedFile> anno_index(gnomad_file);

        //parse gnomad file header
        QByteArrayList gnomad_header;
        foreach(const QByteArray& line, gnomad_file.headers())
        {
            if(line.startsWith("#CHROM_A"))
            {
               gnomad_header = line.trimmed().split('\t');
               break;
            }
        }
        if(gnomad_header.size() <= 0) THROW(FileParseException, "No header line found in gnomad BED file!");

        // Debug
        foreach(QByteArray item, gnomad_header)
        {
            qDebug() << item;
        }

        // get required column indices (-3, because annotation starts at the 4th column
        int i_gnomad_chr2 = gnomad_header.indexOf("CHROM_B") - 3;
        if (i_gnomad_chr2 < 0 ) THROW(FileParseException, "Column 'CHROM_B' not found in gnomAD file!");
        int i_gnomad_start2 = gnomad_header.indexOf("START_B") - 3;
        if (i_gnomad_start2 < 0 ) THROW(FileParseException, "Column 'START_B' not found in gnomAD file!");
        int i_gnomad_end2 = gnomad_header.indexOf("END_B") - 3;
        if (i_gnomad_end2 < 0 ) THROW(FileParseException, "Column 'END_B' not found in gnomAD file!");
        int i_gnomad_type = gnomad_header.indexOf("TYPE") - 3;
        if (i_gnomad_type < 0 ) THROW(FileParseException, "Column 'TYPE' not found in gnomAD file!");

        int i_gnomad_af = gnomad_header.indexOf("AF") - 3;
        if (i_gnomad_af < 0 ) THROW(FileParseException, "Column 'AF' not found in gnomAD file!");
        int i_gnomad_hemi = gnomad_header.indexOf("HEMI") - 3;
        if (i_gnomad_hemi < 0 ) THROW(FileParseException, "Column 'HEMI' not found in gnomAD file!");
        int i_gnomad_hom = gnomad_header.indexOf("HOM") - 3;
        if (i_gnomad_hom < 0 ) THROW(FileParseException, "Column 'HOM' not found in gnomAD file!");

		//af of sub-populations
		int i_gnomad_afr_af = gnomad_header.indexOf("AFR_AF") - 3;
		if (i_gnomad_afr_af < 0 ) THROW(FileParseException, "Column 'AFR_AF' not found in gnomAD file!");
		int i_gnomad_amr_af = gnomad_header.indexOf("AMR_AF") - 3;
		if (i_gnomad_amr_af < 0 ) THROW(FileParseException, "Column 'AMR_AF' not found in gnomAD file!");
		int i_gnomad_eas_af = gnomad_header.indexOf("EAS_AF") - 3;
		if (i_gnomad_eas_af < 0 ) THROW(FileParseException, "Column 'EAS_AF' not found in gnomAD file!");
		int i_gnomad_eur_af = gnomad_header.indexOf("EUR_AF") - 3;
		if (i_gnomad_eur_af < 0 ) THROW(FileParseException, "Column 'EUR_AF' not found in gnomAD file!");
		int i_gnomad_oth_af = gnomad_header.indexOf("OTH_AF") - 3;
		if (i_gnomad_oth_af < 0 ) THROW(FileParseException, "Column 'OTH_AF' not found in gnomAD file!");




		//process BEDPE file
		BedpeFile bedpe_file;
		bedpe_file.load(in);

		// check if annotation already exisits:
		int i_gnomad = bedpe_file.annotationIndexByName("gnomAD", false);
		int i_gnomad_hom_hemi = bedpe_file.annotationIndexByName("gnomAD_hom_hemi", false);
		int i_gnomad_sub = bedpe_file.annotationIndexByName("gnomAD_sub", false);



		// create output buffer and copy comments and header
		QByteArrayList output_buffer;
		output_buffer.append(bedpe_file.comments());
		// get header
		QByteArrayList header = bedpe_file.annotationHeaders();
		// modify header if gene columns not already present
		if (i_gnomad < 0) header.append("gnomAD");
		if (i_gnomad_hom_hemi < 0) header.append("gnomAD_hom_hemi");
		if (i_gnomad_sub < 0) header.append("gnomAD_sub");
		// copy header
		output_buffer << "#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + header.join("\t");


        // debug
        int n_matches = 0;

		for(int i=0; i<bedpe_file.count(); ++i)
		{
            BedpeLine sv = bedpe_file[i];

			//scores
			QByteArrayList gnomad_af;
			QByteArrayList gnomad_hemi;
			QByteArrayList gnomad_hom;
			QByteArrayList gnomad_sub;

			// apply offset
			int sv_s1 = std::max(0, sv.start1() - offset);
			int sv_e1 = sv.end1() + offset;
			int sv_s2 = std::max(0, sv.start2() - offset);
			int sv_e2 = sv.end2() + offset;

			// find matching SV(s)			
			QVector<int> indices;
            if (sv.type() == StructuralVariantType::BND)
			{
				indices = anno_index.matchingIndices(sv.chr1(), sv_s1, sv_e1);
			}
			else
			{
				indices = anno_index.matchingIndices(sv.chr1(), sv_s1, sv_e2);
			}


			foreach (int idx, indices)
			{
                const BedLine& gnomad_line = gnomad_file[idx];

				// filter by type
                if (sv.type() != BedpeFile::stringToType(gnomad_line.annotations().at(i_gnomad_type))) continue;

				// exact match
                if (sv.type() == StructuralVariantType::BND)
                {
                    //check first coos
					if(!gnomad_line.overlapsWith(sv.chr1(), sv_s1, sv_e1)) continue;

                    //check second coos
                    BedLine second_coos = BedLine(Chromosome(gnomad_line.annotations().at(i_gnomad_chr2)),
                                                  Helper::toInt(gnomad_line.annotations().at(i_gnomad_start2)),
                                                  Helper::toInt(gnomad_line.annotations().at(i_gnomad_end2)));
					if(!second_coos.overlapsWith(sv.chr2(), sv_s2, sv_e2)) continue;
                }
                else
                {
                    //check start
					if (sv_s1 > gnomad_line.start() || sv_e1 < gnomad_line.start()) continue;

                    //check end
					if (sv_s2 > gnomad_line.end() || sv_e2 < gnomad_line.end()) continue;
                }

                //debug
//                qDebug() << gnomad_line.toString(true) << gnomad_line.annotations().join("\t");
//                qDebug() << sv.toString();
                n_matches++;

                // extract gnomAD scores
				gnomad_af << gnomad_line.annotations().at(i_gnomad_af);
				gnomad_hom << gnomad_line.annotations().at(i_gnomad_hom);
				gnomad_hemi << gnomad_line.annotations().at(i_gnomad_hemi);

				// parse sub-populations
				QByteArrayList tmp;
				tmp << QByteArray::number(Helper::toDouble(gnomad_line.annotations().at(i_gnomad_afr_af)), 'f', 4);
				tmp << QByteArray::number(Helper::toDouble(gnomad_line.annotations().at(i_gnomad_amr_af)), 'f', 4);
				tmp << QByteArray::number(Helper::toDouble(gnomad_line.annotations().at(i_gnomad_eas_af)), 'f', 4);
				tmp << QByteArray::number(Helper::toDouble(gnomad_line.annotations().at(i_gnomad_eur_af)), 'f', 4);
				tmp << QByteArray::number(Helper::toDouble(gnomad_line.annotations().at(i_gnomad_oth_af)), 'f', 4);

				gnomad_sub << tmp.join(',');
			}


			QByteArray gnomad_af_string;
			QByteArray gnomad_hom_hemi_string;
			QByteArray gnomad_sub_string;

			// check if match found
			if (gnomad_af.size() > 0)
			{
				//warn if multiple matches found
				if (gnomad_af.size() > 1)
				{
					qDebug() << "Multiple matches found for SV " << sv.toString();
				}

				gnomad_af_string = QByteArray::number(Helper::toDouble(gnomad_af.at(0)), 'f', 4);
				gnomad_hom_hemi_string = gnomad_hom.at(0) + "," + gnomad_hemi.at(0);
				gnomad_sub_string = gnomad_sub.at(0);
			}


            QList<QByteArray> annotations = sv.annotations();
			if (i_gnomad > -1)
			{
				annotations[i_gnomad] = gnomad_af_string;
			}
			else
			{
				annotations.append(gnomad_af_string);
			}
			if (i_gnomad_hom_hemi > -1)
			{
				annotations[i_gnomad_hom_hemi] = gnomad_hom_hemi_string;
			}
			else
			{
				annotations.append(gnomad_hom_hemi_string);
			}
			if (i_gnomad_sub > -1)
			{
				annotations[i_gnomad_sub] = gnomad_sub_string;
			}
			else
			{
				annotations.append(gnomad_sub_string);
			}
            sv.setAnnotations(annotations);

			//add annotated line to buffer
            output_buffer << sv.toTsv();
		}


		// open output file and write annotated SVs to file
		QSharedPointer<QFile> cnv_output_file = Helper::openFileForWriting(out, true);
		QTextStream output_stream(cnv_output_file.data());

		foreach (QByteArray line, output_buffer)
		{
			output_stream << line << "\n";
		}
		output_stream.flush();
		cnv_output_file->close();


        // stats
        qDebug() << "Found matches: " << n_matches;
	}

};

#include "main.moc"

int main(int argc, char *argv[])
{
	ConcreteTool tool(argc, argv);
	return tool.execute();
}
