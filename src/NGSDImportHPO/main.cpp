#include "ToolBase.h"
#include "NGSD.h"
#include "Exceptions.h"
#include "Helper.h"
#include "OntologyTermCollection.h"

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
		setDescription("Imports HPO terms and gene-phenotype relations into the NGSD.");
		addInfile("obo", "HPO ontology file from 'http://purl.obolibrary.org/obo/hp.obo'.", false);
		addInfile("anno", "HPO annotations file from 'https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/phenotype_to_genes.txt'", false);

		//optional
		addInfile("omim", "OMIM 'morbidmap.txt' file for additional disease-gene information, from 'https://omim.org/downloads/'.", true);
		addInfile("clinvar", "ClinVar VCF file for additional disease-gene information. Download and unzip from 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2021/clinvar_20210424.vcf.gz' for GRCH37 or 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2021/clinvar_20211212.vcf.gz' for GRCh38.", true);
		addInfile("hgmd", "HGMD phenobase file (Manually download and unzip 'hgmd_phenbase-2021.3.dump').", true);

		// optional (for evidence information):
		addInfile("hpophen", "HPO 'phenotype.hpoa' file for additional phenotype-disease evidence information. Download from https://hpo.jax.org/app/download/annotation", true);
		addInfile("gencc", "gencc 'gencc-submissions.csv' file for additional disease-gene evidence information. Download from https://search.thegencc.org/download", true);
		addInfile("decipher", "G2P 'DDG2P.csv' file for additional gene-disease-phenotype evidence information. Download from https://www.deciphergenomics.org/about/downloads/data", true);

		addFlag("test", "Uses the test database instead of on the production database.");
		addFlag("force", "If set, overwrites old data.");
		addFlag("debug", "Enables debug output");

		changeLog(2021,12,22, "Added support for GenCC and DECIPHER.");
		changeLog(2020, 7, 7, "Added support of HGMD gene-phenotype relations.");
		changeLog(2020, 3, 5, "Added support for new HPO annotation file.");
		changeLog(2020, 3, 9, "Added optimization for hpo-gene relations.");
		changeLog(2020, 3, 10, "Removed support for old HPO annotation file.");
		changeLog(2020, 7, 6, "Added support for HGMD phenobase file.");
	}

	/// simple sruct to keep a set of source databases
	struct SourceDetails {
		QList<PhenotypeSource::Source> sources;
		QStringList original_evidence;
		QList<PhenotypeEvidence::Evidence> translated_evidence;

		SourceDetails()
		{
			sources = QList<PhenotypeSource::Source>();
		}

		SourceDetails(const QByteArray& s, const QByteArray& original_evi, PhenotypeEvidence::Evidence translated_evi=PhenotypeEvidence::Evidence::NA)
		{
			sources = QList<PhenotypeSource::Source>();
			original_evidence = QStringList();
			translated_evidence = QList<PhenotypeEvidence::Evidence>();

			sources.append(PhenotypeSource::sourceFromString(s));
			original_evidence.append(QString(original_evi));
			translated_evidence.append(translated_evi);
		}

		bool contains(const PhenotypeSource::Source& s)
		{
			return sources.contains(s);
		}

		bool contains(const QByteArray& s)
		{
			return sources.contains(PhenotypeSource::sourceFromString(s));
		}

		bool contains(const QString& s)
		{
			return sources.contains(PhenotypeSource::sourceFromString(s));
		}

		int getIndexOfSource(const PhenotypeSource::Source& s)
		{
			return sources.indexOf(s);
		}

		void append(const PhenotypeSource::Source& s,const QString& original_evi, PhenotypeEvidence::Evidence translated_evi)
		{
			sources.append(s);
			original_evidence.append(original_evi);
			translated_evidence.append(translated_evi);
		}

		void append(const QByteArray& s, const QString& original_evi, PhenotypeEvidence::Evidence translated_evi)
		{
			sources.append(PhenotypeSource::sourceFromString(s));
			original_evidence.append(original_evi);
			translated_evidence.append(translated_evi);
		}

		void unite(const SourceDetails& second)
		{
			for (int i=0; i<second.sources.count(); i++)
			{
				if (contains(second.sources[i]))
				{
					int idx = getIndexOfSource(second.sources[i]);
					if ((int) second.translated_evidence[i] > (int) translated_evidence[idx])
					{
						original_evidence[idx] = second.original_evidence[i];
						translated_evidence[idx] = second.translated_evidence[i];
					}
				}
				else
				{
					append(second.sources[i], second.original_evidence[i], second.translated_evidence[i]);
				}
			}
		}

		QString toCsvString() const
		{
			QString s = "";

			for (int i=0; i<sources.count(); i++)
			{
				if (i == 0)
				{
					s= "(" + PhenotypeSource::sourceToString(sources[i]) + ", " + original_evidence[i] + ", " + PhenotypeEvidence::evidenceToString(translated_evidence[i]) + ")";
				}
				else
				{
					s += "; (" + PhenotypeSource::sourceToString(sources[i]) + ", " + original_evidence[i] + ", " + PhenotypeEvidence::evidenceToString(translated_evidence[i]) + ")";
				}
			}
			return s;
		}
	};

	// struct used to find the first instance of where a relation was parsed from. Only used for debugging and testing
	struct ExactSources
	{
		QString disease2gene;
		QString term2disease;
		QString term2gene;

		ExactSources():
		disease2gene("")
		, term2disease("")
		, term2gene("")
		{
		}

		ExactSources(const QString& disease2phenotype, const QString& phenotype2gene, const QString& term2gene):
		disease2gene(disease2phenotype)
		, term2disease(phenotype2gene)
		, term2gene(term2gene)
		{
		}

		void combine(const ExactSources& src)
		{
			if (disease2gene.length() == 0)
			{
				disease2gene = src.disease2gene;
			}

			if (term2disease.length() == 0)
			{
				term2disease = src.term2disease;
			}

			if (term2gene.length() == 0)
			{
				term2gene = src.term2gene;
			}
		}

		QString toString() const
		{
			return QString("Exact Sources: disease2gene \t  '%1', \t term2disease \t '%2', \t term2gene \t '%3'").arg(disease2gene, term2disease, term2gene);
		}
	};

	// struct to annotate a gene/disease/phenotype with a source db and an evidence strength
	struct AnnotatedItem
	{
		QByteArray item;
		SourceDetails src;
		PhenotypeEvidence::Evidence evi;
		ExactSources exactSources; // for debuging and testing

		AnnotatedItem()
		{
		}

		AnnotatedItem(const QByteArray& item, const QByteArray& s, const QByteArray& original_evi, PhenotypeEvidence::Evidence evi):
			item(item)
		  , evi(evi)
		{
			src = SourceDetails(s, original_evi, evi);
		}

		AnnotatedItem(const QByteArray& item, SourceDetails src, PhenotypeEvidence::Evidence evi):
			item(item)
		  , src(src)
		  , evi(evi)
		{
		}

		AnnotatedItem(const QByteArray& item, SourceDetails src, PhenotypeEvidence::Evidence evi, ExactSources exactSources):
			item(item)
		  , src(src)
		  , evi(evi)
		  , exactSources(exactSources)
		{
		}

		bool operator==(const AnnotatedItem& other)
		{
			return item == other.item;
		}
	};

	/// a QList of AnnotatedItems with some convenience insertion methods
	class AnnotatedList
	{
		public:
			void add(const QByteArray& item, const QByteArray& source, const QByteArray& original_evi, PhenotypeEvidence::Evidence evidence=PhenotypeEvidence::NA, ExactSources exactSource=ExactSources())
			{
				add(item, SourceDetails(source, original_evi, evidence), evidence, exactSource);
			}

			void add(const QByteArray& item, SourceDetails source, PhenotypeEvidence::Evidence evidence=PhenotypeEvidence::NA, ExactSources exactSource=ExactSources())
			{
				if (hash.contains(item))
				{
					AnnotatedItem& present_item = hash[item];
					present_item.src.unite(source);
					if ((int) present_item.evi < (int) evidence)
					{
						present_item.evi = evidence;
						present_item.exactSources.combine(exactSource);
					}

				}
				else
				{
					hash.insert(item, AnnotatedItem(item, source, evidence, exactSource));
				}
			}

			QList<AnnotatedItem> items() const
			{
				return hash.values();
			}

		protected:
			QHash<QByteArray, AnnotatedItem> hash;

	};

	QHash<QByteArray, int> importHpoOntology(const NGSD& db)
	{
		QTextStream out(stdout);
		//prepare SQL queries
		SqlQuery qi_term = db.getQuery();
		qi_term.prepare("INSERT INTO hpo_term (hpo_id, name, definition, synonyms) VALUES (:0, :1, :2, :3);");
		SqlQuery qi_parent = db.getQuery();
		qi_parent.prepare("INSERT INTO hpo_parent (parent, child) VALUES (:0, :1);");

		QHash<QByteArray, int> id2ngsd;
		OntologyTermCollection terms(getInfile("obo"), true);
		for (int i=0; i<terms.size(); ++i)
		{
			const OntologyTerm& term = terms.get(i);

			qi_term.bindValue(0, term.id());
			qi_term.bindValue(1, term.name());
			qi_term.bindValue(2, term.definition());
			qi_term.bindValue(3, term.synonyms().count()==0 ? "" : term.synonyms().join('\n'));
			qi_term.exec();

			id2ngsd.insert(term.id(), qi_term.lastInsertId().toInt());
		}
		out << "Imported " << db.getValue("SELECT COUNT(*) FROM hpo_term").toInt() << " non-obsolete HPO terms." << endl;

		//insert parent-child relations between (valid) terms
		for (int i=0; i<terms.size(); ++i)
		{
			const OntologyTerm& term = terms.get(i);

			int c_db = id2ngsd.value(term.id(), -1);
			if (c_db==-1) continue;

			foreach(const QByteArray& p_id, term.parentIDs())
			{
				int p_db = id2ngsd.value(p_id, -1);
				if (p_db==-1) continue;

				qi_parent.bindValue(0, p_db);
				qi_parent.bindValue(1, c_db);
				qi_parent.exec();
			}
		}
		out << "Imported " << db.getValue("SELECT COUNT(*) FROM hpo_parent").toInt() << " parent-child relations between terms." << endl;

		return id2ngsd;
	}

	void parseHpoPhen(QHash<QByteArray, int> id2ngsd, QHash<int, AnnotatedList>& term2diseases)
	{
		if (getInfile("hpophen") == "") return;

		int lineCount = 0;
		int added = 0;

		QTextStream out(stdout);
		if (getFlag("debug")) out << "Starting analysis of hpophen file\n";

		// parse phenotype.hpoa file for evidence information
		QSharedPointer<QFile> fp = Helper::openFileForReading(getInfile("hpophen"));

		while(! fp->atEnd())
		{
			lineCount++;
			QByteArray line = fp->readLine();

			if (line.startsWith('#')) continue;

			QByteArrayList parts = line.split('\t');

			if (parts[2].length() > 0) continue; // Qualifier Not: term NOT associated to the disease

			QByteArray term = parts[3].trimmed();
			QByteArray disease = parts[4].trimmed();
			QByteArray evidence = parts[5].trimmed();

			int term_id = id2ngsd.value(term, -1);

			if (term_id == -1)
			{
				if (getFlag("debug"))
				{
					out << "Term not found in id2ngsd: " << term << "\n";
				}
				continue;
			}
			else
			{
				ExactSources e_src = ExactSources();
				e_src.term2disease = QString("hpoPhen line ") + QString::number(lineCount);
				term2diseases[term_id].add(disease, "HPO", evidence,  PhenotypeEvidence::translateHpoEvidence(evidence), e_src);
				added++;
				if (getFlag("debug"))
				{
					out << "Imported term2disease relation:\t" << term << "-" << disease << ":\t" << evidence << "\t fin_evi:\t" << PhenotypeEvidence::evidenceToString(PhenotypeEvidence::translateHpoEvidence(evidence)) << "\n";
				}
			}
		}
		out << "Imported " << added << " term-disease relations.\n";
		fp->close();
	}

	void parseDecipher(NGSD& db, const QHash<QByteArray, int>& id2ngsd, QHash<QByteArray, AnnotatedList>& disease2genes, QHash<int, AnnotatedList>& term2diseases, QHash<int, AnnotatedList>& term2genes)
	{
		if (getInfile("decipher") == "") return;
		bool debug = getFlag("debug");
		QTextStream out(stdout);
		if (debug) out << "Parsing Decipher...\n";
		int countT2D = 0;
		int countD2G = 0;
		int countT2G = 0;

		QSharedPointer<QFile> fp = Helper::openFileForReading(getInfile("decipher"));

		QByteArray line = fp->readLine();
		//"gene symbol","gene mim","disease name","disease mim","confidence category","allelic requirement","mutation consequence",phenotypes,"organ specificity list",pmids,panel,"prev symbols","hgnc id","gene disease pair entry date","cross cutting modifier","mutation consequence flag"

		QList<QByteArray> non_hgc_genes;
		QSet<QByteArray> bad_hpo_terms;
		QByteArray source = "Decipher";
		int lineCount = 0;
		QRegExp mim_exp("([0-9]{6})");
		while(! fp->atEnd())
		{
			lineCount++;
			line = fp->readLine().trimmed();
			QByteArrayList parts = line.split(',');

			// merge parts of strings that contained commas:

			parts = reconstructStrings(parts);

			QByteArray gene = parts[0].trimmed();
			QByteArray disease_num = parts[3].trimmed();
			QByteArray disease = "OMIM:" + parts[3].trimmed();
			QByteArray decipher_evi = parts[4].trimmed();
			PhenotypeEvidence::Evidence evidence = PhenotypeEvidence::translateDecipherEvidence(decipher_evi);
			QByteArrayList hpo_terms = parts[7].trimmed().split(';');

			//verify information
			int gene_db_id = db.geneToApprovedID(gene);

			if (gene_db_id == -1)
			{
				non_hgc_genes.append(gene);
				// add term2disease relations
				foreach (const QByteArray& term, hpo_terms)
				{
					int term_db_id = id2ngsd.value(term, -1);
					if (term_db_id != -1 && mim_exp.indexIn(disease_num) != -1)
					{
						ExactSources e_src = ExactSources();
						e_src.term2disease = QString("Decipher line") + QString::number(lineCount);
						term2diseases[term_db_id].add(disease, source, decipher_evi,  evidence, e_src);
						countT2D++;
						if (debug) out << "Deciper\tTERM2DISEASE\tTERM,DISEASE,GENE\t" << term << "\t" << disease << "\t''\t" << "\tSource:\t" << e_src.term2disease <<"\n";
					}
					else
					{
						bad_hpo_terms.insert(term);
					}
				}
			}
			else
			{
				QByteArray approved_gene_symbol = db.geneSymbol(gene_db_id);

				foreach (const QByteArray& term, hpo_terms)
				{
					int term_db_id = id2ngsd.value(term, -1);
					if (term_db_id != -1)
					{
						ExactSources e_src = ExactSources();
						e_src.term2gene = QString("Decipher line") + QString::number(lineCount);
						term2genes[term_db_id].add(approved_gene_symbol, source, decipher_evi, evidence, e_src);
						countT2G++;
						if (debug) out << "Deciper\tTERM2GENE\tTERM,DISEASE,GENE\t" << term << "\t''\t" << gene << "\tSource:\t" << e_src.term2gene << "\tapproved_gene_symbol:\t" << approved_gene_symbol << "\n";

						if (mim_exp.indexIn(disease_num) != -1)
						{
							e_src = ExactSources();
							e_src.term2disease = QString("Decipher line") + QString::number(lineCount);
							term2diseases[term_db_id].add(disease, source, decipher_evi, evidence, e_src);
							countT2D++;
							if (debug) out << "Deciper\tTERM2DISEASE\tTERM,DISEASE,GENE\t" << term << "\t" << disease << "\t''\t" << "\tSource:\t" << e_src.term2disease <<"\n";
						}

					}
					else
					{
						bad_hpo_terms.insert(term);
					}
				}
				if (mim_exp.indexIn(disease_num) != -1)
				{
					ExactSources e_src = ExactSources();
					e_src.term2disease = QString("Decipher line") + QString::number(lineCount);
					disease2genes[disease].add(approved_gene_symbol, source, decipher_evi, evidence, e_src);
					countD2G++;
					if (debug) out << "Deciper\tTERM2GENE\tTERM,DISEASE,GENE\t" << "''\t" << disease << "\t" << gene << "\tSource:\t" << e_src.term2gene << "\tapproved_gene_symbol:\t" << approved_gene_symbol << "\n";
				}
			}
		}
		fp->close();
		out << "Finished Deciper parsing. Imported:\n" << countD2G << "\tDisease to gene relations\n" << countT2D << "\tHPO term to disease relations\n" << countT2G << "\tHPO term to gene relations.\n";

	}

	void parseGenCC(NGSD& db, QHash<QByteArray, AnnotatedList>& disease2genes)
	{
		if (getInfile("gencc") == "") return;

		// parse gencc_submission.csv file for evidence information
		QSharedPointer<QFile> fp = Helper::openFileForReading(getInfile("gencc"));
		QByteArray line = fp->readLine(); // header
		QByteArray source = "GenCC";
		int count =0;
		int lineCount =0;
		while(! fp->atEnd())
		{
			lineCount++;
			line = fp->readLine().trimmed();
			while ( ! line.endsWith('"')) //some strings contain newlines..
			{
				lineCount++;
				line.append(fp->readLine().trimmed());
			}

			QByteArrayList parts = line.split(',');

			parts = reconstructStrings(parts, 30);

			QByteArray gene_symbol = parts[2].replace('"', ' ').trimmed();
			QByteArray disease = parts[5].replace('"', ' ').trimmed(); // OMIM:XXXXXX, MONDO:XXXXXXX, Orphanet:XXXXX needs mapping from Orphanet and Mondo to Omim
			QByteArray gencc_evi = parts[8].replace('"', ' ').trimmed();
			PhenotypeEvidence::Evidence evidence = PhenotypeEvidence::translateGenccEvidence(gencc_evi);

			if (evidence == PhenotypeEvidence::NA || evidence == PhenotypeEvidence::AGAINST)
			{
				continue;
			}

			if ( ! disease.startsWith("OMIM"))
			{
				continue;
			}

			int gene_db_id = db.geneToApprovedID(gene_symbol);
			if (gene_db_id == -1) continue;
			ExactSources e_src = ExactSources();
			e_src.disease2gene = QString("GenCC line") + QString::number(lineCount);
			disease2genes[disease].add(db.geneSymbol(gene_db_id), source, gencc_evi, evidence, e_src);
			count++;
		}
		fp->close();

		QTextStream out(stdout);
		out << "Imported " << count << " disease gene relations from GenCC" << endl;
	}

	QByteArrayList reconstructStrings(const QByteArrayList& parts, int expected_size=-1)
	{
		// if parts size bigger than expected, try to reconstruct strings that were split:
		if (parts.length() > expected_size)
		{
			QByteArrayList cleaned_parts = QByteArrayList();

			for (int i=0; i<parts.length(); i++)
			{
				if (parts[i].startsWith('"') && ( ! parts[i].endsWith('"'))) // starts with " but doesn't end with "
				{
					QByteArray combined_part = parts[i];
					do
					{
						i++;
						combined_part.append(parts[i]);

					}
					while (! parts[i].endsWith('"'));

					cleaned_parts.append(combined_part);
				}
				else
				{
					cleaned_parts.append(parts[i]);
				}

			}
			return cleaned_parts;
		}

		return parts;
	}

	virtual void main()
	{
		//init
		NGSD db(getFlag("test"));
		QTextStream out(stdout);
		bool debug = getFlag("debug");

		//check tables exist
		db.tableExists("hpo_term");
		db.tableExists("hpo_parent");
		db.tableExists("hpo_genes");

		//check if gene table exists and contains HGNC genes:
		db.tableExists("gene");
		SqlQuery test_gene_table = db.getQuery();
		test_gene_table.exec("SELECT count(*) FROM gene;");
		while (test_gene_table.next())
		{
			if (test_gene_table.value(0) == 0)
			{
				THROW(DatabaseException, "Table 'gene' is empty. Please import HGNC database before importing HPO.")
			}
		}

		//clear tables if not empty
		if (!db.tableEmpty("hpo_term") || !db.tableEmpty("hpo_parent") || !db.tableEmpty("hpo_genes"))
		{
			if (getFlag("force"))
			{
				db.clearTable("hpo_genes");
				db.clearTable("hpo_parent");
				db.clearTable("hpo_term");
			}
			else
			{
				THROW(DatabaseException, "Tables already contain data! Use '-force' to overwrite old data!");
			}
		}

		// import HPO Ontology to DB
		QHash<QByteArray, int> id2ngsd = importHpoOntology(db);

		QHash<int, AnnotatedList> term2genes;
		QHash<int, AnnotatedList > term2diseases;
		QHash<QByteArray, AnnotatedList> disease2genes;

		// parse Evidence files if provided
		// parse g2pDDG2P_11_11_2021.csv file
		parseDecipher(db, id2ngsd, disease2genes, term2diseases, term2genes);
		// parse gencc-submissions.csv file
		parseGenCC(db, disease2genes);
		// parse phenotype.hpoa file
		parseHpoPhen(id2ngsd, term2diseases);


		//parse term-disease and disease-gene relations from HPO
		QSharedPointer<QFile> fp = Helper::openFileForReading(getInfile("anno"));
		QSet<QByteArray> non_hgnc_genes;
		PhenotypeList inheritance_terms = db.phenotypeChildTerms(db.phenotypeIdByAccession("HP:0000005"), true); //Mode of inheritance
		int lineCount = 0;
		QString exactSource;
		while(!fp->atEnd())
		{
			lineCount++;
			exactSource = QString("Anno line ") + QString::number(lineCount);
			QByteArray line =  fp->readLine();
			QByteArrayList parts =line.split('\t');

			if (parts.count()<7) continue;

			// parse line
			QByteArray disease = parts[6].trimmed();
			QByteArray gene = parts[3].trimmed();
			QByteArray term_accession = parts[0].trimmed();

			int gene_db_id = db.geneToApprovedID(gene);
			int term_db_id = id2ngsd.value(term_accession, -1);

			if (term_db_id!=-1)
			{
				if (inheritance_terms.containsAccession(term_accession))
				{
					if (gene_db_id!=-1)
					{
						if (debug) out << "HPO-GENE: " << term_accession << " - " << gene << "\n";

						ExactSources e_src = ExactSources();
						e_src.term2gene = exactSource;
						term2genes[term_db_id].add(db.geneSymbol(gene_db_id), "HPO", "", PhenotypeEvidence::NA, e_src);
					}
				}
				else
				{
					if (debug) out << "HPO-DISEASE: " << term_accession << " - " << disease << "\n";

					ExactSources e_src = ExactSources();
					e_src.term2disease = exactSource;
					term2diseases[term_db_id].add(disease, "HPO", "", PhenotypeEvidence::NA, e_src);
				}
			}

			if (gene_db_id!=-1)
			{
				if (debug) out << "DISEASE-GENE (HPO): " << disease << " - " << db.geneSymbol(gene_db_id) << "\n";

				ExactSources e_src = ExactSources();
				e_src.disease2gene = exactSource;
				disease2genes[disease].add(db.geneSymbol(gene_db_id), "HPO", "", PhenotypeEvidence::NA, e_src);
			}
			else
			{
				non_hgnc_genes << gene;
			}
		}
		fp->close();

		foreach(const QByteArray& gene, non_hgnc_genes)
		{
			out << "Skipped gene '" << gene << "' because it is not an approved HGNC symbol!" << endl;
		}

		//parse disease-gene relations from OMIM
		QString omim_file = getInfile("omim");
		lineCount = 0;
		int count = 0;
		if (omim_file!="")
		{
			if (debug) out << "Parsing OMIM file...\n";
			//parse disease-gene relations
			int c_skipped_invalid_gene = 0;
			fp = Helper::openFileForReading(omim_file);
			QRegExp mim_exp("([0-9]{6})");
			QRegExp evi_exp("(\\([1-4]{1}\\))");

			while(!fp->atEnd())
			{
				lineCount++;
				QByteArrayList parts = fp->readLine().trimmed().split('\t');
				if (parts.count()<4) continue;

				QByteArray pheno = parts[0].trimmed();
				QByteArrayList genes = parts[1].split(',');
				QByteArray mim_number = parts[2].trimmed(); // mim number for gene
				QByteArray omim_evi = "";

				if (mim_exp.indexIn(pheno)!=-1)
				{
					mim_number = mim_exp.cap().toLatin1(); // mim number for phenotype
				}

				if (evi_exp.indexIn(pheno) != -1)
				{
					omim_evi = evi_exp.cap().toLatin1(); // evidence for relation
				}

				foreach(QByteArray gene, genes)
				{
					//make sure the gene symbol is approved by HGNC
					gene = gene.trimmed();
					int approved_id = db.geneToApprovedID(gene);
					if (approved_id==-1)
					{
						if (debug) out << "Skipped gene '" << gene << "' because it is not an approved HGNC symbol!\n";
						++c_skipped_invalid_gene;
						continue;
					}

					if (debug) out << "DISEASE-GENE (OMIM): OMIM:" << mim_number << " - " << db.geneSymbol(approved_id) << "\n";

					ExactSources e_src = ExactSources();
					e_src.disease2gene = QString("OMIM line ") + QString::number(lineCount);
					disease2genes["OMIM:"+mim_number].add(db.geneSymbol(approved_id), "OMIM", omim_evi, PhenotypeEvidence::translateOmimEvidence(omim_evi), e_src);
					count++;
				}
			}
			fp->close();
			out << "Imported " << count << " disease to genes relations.\n";
		}

		//parse disease-gene relations from ClinVar
		QString clinvar_file = getInfile("clinvar");
		if (clinvar_file!="")
		{
			if (debug) out << "Prasing ClinVar..." << endl;
			//parse disease-gene relations
			int c_skipped_invalid_gene = 0;
			fp = Helper::openFileForReading(clinvar_file);
			lineCount = 0;
			while(!fp->atEnd())
			{
				lineCount++;
				QByteArray line = fp->readLine().trimmed();
				if (!line.contains("CLNSIG=Pathogenic") && !line.contains("CLNSIG=Likely_pathogenic")) continue;

				QByteArrayList parts = line.split('\t');
				if (parts.count()<8) continue;

				//parse gene/disease info from from INFO field
				GeneSet genes;
				QByteArrayList diseases;
				QByteArrayList hpos;
				QByteArrayList info_parts = parts[7].split(';');
				foreach(const QByteArray& part, info_parts)
				{
					if (part.startsWith("GENEINFO="))
					{
						QByteArrayList geneinfo_parts = part.mid(9).split('|');
						foreach(QByteArray geneinfo, geneinfo_parts)
						{
							int colon_idx = geneinfo.indexOf(':');
							if (colon_idx==-1) continue;
							genes << geneinfo.left(colon_idx);
						}
					}
					if (part.startsWith("CLNDISDB=") || part.startsWith("CLNDISDBINCL="))
					{
						QByteArrayList disease_parts = part.mid(part.indexOf('=')+1).replace(',', '|').split('|');
						foreach(QByteArray disease_part, disease_parts)
						{
							if (disease_part.startsWith("OMIM:"))
							{
								diseases << disease_part.trimmed();
							}
							else if (disease_part.startsWith("Orphanet:ORPHA"))
							{
								diseases << disease_part.replace("Orphanet:ORPHA", "ORPHA:").trimmed();
							}
							else if (disease_part.startsWith("Human_Phenotype_Ontology:"))
							{
								hpos << disease_part.replace("Human_Phenotype_Ontology:", "").trimmed();
							}
						}
					}
				}
				if (genes.isEmpty() || (diseases.isEmpty() && hpos.isEmpty())) continue;

				foreach(const QByteArray& gene, genes)
				{
					//make sure the gene symbol is approved by HGNC
					int approved_id = db.geneToApprovedID(gene);
					if (approved_id==-1)
					{
						if (debug) out << "Skipped gene '" << gene << "' because it is not an approved HGNC symbol!" << endl;
						++c_skipped_invalid_gene;
						continue;
					}
					QByteArray gene_approved = db.geneSymbol(approved_id);

					foreach(const QByteArray& disease, diseases)
					{
						if (debug) out << "DISEASE-GENE (ClinVar): " << disease << " - " << gene_approved << "\n";

						ExactSources e_src = ExactSources();
						e_src.disease2gene = QString("ClinVar line ") + QString::number(lineCount);
						disease2genes[disease].add(gene_approved, "ClinVar", "", PhenotypeEvidence::NA, e_src);
					}
					foreach(const QByteArray& hpo, hpos)
					{
						if (debug) out << "HPO-GENE (ClinVar): " << hpo << " - " << gene_approved << "\n";
						int term_db_id = id2ngsd.value(hpo, -1);
						if (term_db_id != -1)
						{
							ExactSources e_src = ExactSources();
							e_src.term2gene = QString("ClinVar line ") + QString::number(lineCount);
							term2genes[term_db_id].add(gene_approved, "ClinVar", "", PhenotypeEvidence::NA, e_src);
						}
					}
				}
			}
			fp->close();
		}

		// parse hpo-gene relations from HGMD (Phenobase dbdump file):
		QString hgmd_file = getInfile("hgmd");
		if(hgmd_file != "")
		{
			if (debug) out << "Parsing HGMD Phenobase dump file...\n" << endl;
			// define look-up tables
			QMultiMap<int, QByteArray> phenid2gene_mapping = QMap<int, QByteArray>();
			QMultiMap<QByteArray,int> cui2phenid_mapping = QMap<QByteArray,int>();
			QMultiMap<QByteArray,QByteArray> hpo2cui_mapping = QMap<QByteArray,QByteArray>();

			QSharedPointer<QFile> fp = Helper::openFileForReading(hgmd_file);
			int line_number = 0;
			while(!fp->atEnd())
			{
				line_number++;
				// show progress
				if(debug && line_number%100 == 0) out << "\tparsed " << line_number << " lines..." << endl;
				QByteArray line = fp->readLine().trimmed();
				if (line.isEmpty()) continue;

				// parse concept table
				if(line.startsWith("INSERT INTO `concept` VALUES "))
				{
				   // parse insert line
				   QString value_string = line.mid(31);
				   value_string.chop(3);
				   QStringList tuples = value_string.split("'),('");

				   // parse tuples
				   foreach(const QString& tuple, tuples)
				   {
					   QStringList tuple_entries = tuple.split("','");

					   // check if tuple size is correct
					   if(tuple_entries.size() != 10)
					   {
						   THROW(FileParseException, "Invalid number of columns in INSERT Statement in line " + QString::number(line_number) + "\n" + tuple);
					   }

					   // ignore all non HPO entries:
					   if(tuple_entries.at(2).trimmed() != "HPO") continue;

					   QByteArray cui = tuple_entries.at(0).toUtf8();
					   QByteArray hpo = tuple_entries.at(3).toUtf8();

					   // skip already stored realations (due to old/different description/synonyms)
					   if(hpo2cui_mapping.contains(hpo, cui)) continue;

					   hpo2cui_mapping.insert(hpo, cui);
				   }
				}
				else if(line.startsWith("INSERT INTO `hgmd_mutation` VALUES "))
				{
					// parse insert line
					QString value_string = line.mid(36);
					value_string.chop(2);
					QStringList tuples = value_string.split("),(");

					// parse tuples
					foreach(const QString& tuple, tuples)
					{
						QStringList tuple_entries = tuple.split(",");

						// check if tuple size is correct
						if(tuple_entries.size() != 3)
						{
							THROW(FileParseException, "Invalid number of columns in INSERT Statement in line " + QString::number(line_number) + "\n" + tuple);
						}

						//parse gene name
						QByteArray gene_name = tuple_entries.at(1).toUtf8();
						int phen_id = Helper::toInt(tuple_entries.at(2).toUtf8(), "phen_id", QString::number(line_number));
						//remove quotes
						gene_name.remove(0, 1);
						gene_name.remove((gene_name.size() - 1), 1);

						// skip already stored mappings to save memory
						if(phenid2gene_mapping.contains(phen_id, gene_name)) continue;

						phenid2gene_mapping.insert(phen_id, gene_name);
					}
				}
				else if(line.startsWith("INSERT INTO `phenotype_concept` VALUES "))
				{
					// parse insert line
					QString value_string = line.mid(40);
					value_string.chop(2);
					QStringList tuples = value_string.split("),(");

					// parse tuples
					foreach(const QString& tuple, tuples)
					{
						QStringList tuple_entries = tuple.split(",");

						// check if tuple size is correct
						if(tuple_entries.size() != 3)
						{
							THROW(FileParseException, "Invalid number of columns in INSERT Statement in line " + QString::number(line_number) + "\n" + tuple);
						}

						//parse cui
						QByteArray cui = tuple_entries.at(2).toUtf8();
						int phen_id = Helper::toInt(tuple_entries.at(0).toUtf8(), "phen_id", QString::number(line_number));
						//remove quotes
						cui.remove(0, 1);
						cui.remove((cui.size() - 1), 1);

						// skip already stored mappings to save memory
						if(cui2phenid_mapping.contains(cui, phen_id)) continue;

						cui2phenid_mapping.insert(cui, phen_id);
					}
				}
			}


			// create hpo-gene relation
			if (debug) out << "Creating HPO-gene relation from HGMD file..." << endl;

			int hgmd_skipped_invalid_gene = 0;
			int hgmd_genes_added = 0;
			int i = 0;
			foreach (const QByteArray& hpo, hpo2cui_mapping.keys())
			{
				i++;

				// get the corresponding cui for each hpo term:
				foreach(const QByteArray& cui, hpo2cui_mapping.values(hpo))
				{
					// get all phenotype ids for given cui
					if(!cui2phenid_mapping.contains(cui))
					{
					   if (debug) out << "No phenotype id found for CUI '" << cui << "' (HGMD file)!" << endl;
					   continue;
					}
					foreach(int phen_id, cui2phenid_mapping.values(cui))
					{
						// get all genes for a given phenotype id
						if(!phenid2gene_mapping.contains(phen_id))
						{
							if (debug) out << "No gene found for phenotype id " << QByteArray::number(phen_id) << " (HGMD file)!" << endl;
							continue;
						}
						foreach(const QByteArray& gene, phenid2gene_mapping.values(phen_id))
						{
							//make sure the gene symbol is approved by HGNC
							int approved_id = db.geneToApprovedID(gene);
							if (approved_id==-1)
							{
									if (debug) out << "Skipped gene '" << gene << "' because it is not an approved HGNC symbol (HGMD file)!" << endl;
									++hgmd_skipped_invalid_gene;
									continue;
							}
							QByteArray gene_approved = db.geneSymbol(approved_id);

							// add gene to hpo list:
							if (debug) out << "HPO-GENE (HGMD): " << hpo << " - " << gene_approved << endl;
							hgmd_genes_added++;
							int term_db_id = id2ngsd.value(hpo, -1);
							if (term_db_id != -1)
							{
								ExactSources e_src = ExactSources();
								e_src.term2gene = QString("HGMD unknown line");
								term2genes[term_db_id].add(gene_approved, "HGMD", "", PhenotypeEvidence::NA, e_src); // is there some evidence in the file that could be parsed?
							}
						}
					}
				}

				// show progress
				if(debug && (i%1000 == 0))
				{
					out << "\t" << i << " of " << hpo2cui_mapping.keys().size() << "hpo terms parsed \n";
				}
			}
		}

		// import gathered data:
		out << "Gathering all term2gene relations" << endl;

		foreach (int term_id, term2diseases.keys())
		{
			foreach (const AnnotatedItem& disease, term2diseases[term_id].items())
			{
				foreach (const AnnotatedItem& gene, disease2genes[disease.item].items())
				{
					// if one of the evidencess is NA take the other one. If both have a value take the lower ranked one.
					PhenotypeEvidence::Evidence evi;
					if (disease.evi == PhenotypeEvidence::NA)
					{
						evi = gene.evi;
					}
					else if (gene.evi == PhenotypeEvidence::NA)
					{
						evi = disease.evi;
					}
					else
					{
						evi = (int) disease.evi < (int) gene.evi ? disease.evi : gene.evi;
					}

					SourceDetails src = SourceDetails(); // list all the combined sources?
					src.unite(disease.src);
					src.unite(gene.src);
					ExactSources exactSource = disease.exactSources;
					exactSource.combine(gene.exactSources);

					term2genes[term_id].add(gene.item, src, evi, exactSource);
				}
			}
		}

		out << "Starting import into NGSD" << endl;
		// build insert statements with 10000 tuples each (big insertions are way faster than single element insert statements)
		QStringList tuples;
		foreach (int term_id, term2genes.keys())
		{
			foreach (const AnnotatedItem& gene, term2genes[term_id].items())
			{
				if (getFlag("debug"))
				{
					out << "Gene:\t" << gene.item << "\tHPO term id:\t" << term_id << "\t" << "final evidence:\t" << PhenotypeEvidence::evidenceToString(gene.evi) << "  \torigin:\t" << gene.exactSources.toString() << "\n";
				}
				tuples << QString("(%1, '%2', '%3', '%4')").arg(QString::number(term_id), QString(gene.item), gene.src.toCsvString(), PhenotypeEvidence::evidenceToString(gene.evi));
			}
		}
		//import
		for (int i=0; i<tuples.count(); i+=10000)
		{
			db.getQuery().exec("INSERT INTO hpo_genes (hpo_term_id, gene, details, evidence) VALUES " + tuples.mid(i, 10000).join(", ") +";");
		}
		out << "Overall imported term-gene relations: " << db.getValue("SELECT COUNT(*) FROM hpo_genes").toInt() << endl;

		out << "Optimizing term-gene relations...\n";
		out << "(removing genes which are present in leaf nodes from parent node)" << endl;

		Phenotype root = Phenotype("HP:0000001", "All");
		int removed_genes = 0;
		optimizeHpoGeneTable(root, db, id2ngsd, removed_genes);

		// compute import stats

		// get first level of subtrees:
		PhenotypeList subtree_roots = db.phenotypeChildTerms(db.phenotypeIdByAccession(root.accession()), false);
		QList<PhenotypeList> subtrees;
		foreach (const Phenotype& pt, subtree_roots)
		{
			subtrees.append(db.phenotypeChildTerms(db.phenotypeIdByAccession(pt.accession()), true));
		}
		QVector<int> subtree_counts(subtree_roots.count());

		//calulate stats:
		QStringList hpo_terms = db.getValues("SELECT ht.hpo_id, hg.gene FROM hpo_genes hg INNER JOIN hpo_term ht ON hg.hpo_term_id = ht.id");
		foreach (const QString& hpo_term, hpo_terms)
		{
			QByteArray pt = hpo_term.toUtf8();
			for (int i = 0; i < subtree_roots.count(); ++i)
			{
				if (subtrees[i].containsAccession(pt)) subtree_counts[i]++;
			}
		}

		out << "Imported HPO-Gene relations: \n";
		out << " Overall:\t" << hpo_terms.size() << "\n";
		for (int i = 0; i < subtree_roots.count(); ++i)
		{
			out << " " << subtree_roots[i].name() << ":\t" << subtree_counts.at(i) << "\n";
		}

		out << removed_genes << " duplicate genes removed during optimization" << endl;
	}

	void optimizeHpoGeneTable(const Phenotype& root, NGSD& db, const QHash<QByteArray, int>& pt2id, int& removed_genes)
	{
		// get all child nodes
		int root_id = db.phenotypeIdByAccession(root.accession());
		PhenotypeList children = db.phenotypeChildTerms(root_id, false);

		// abort if leaf node
		if (children.count() == 0) return;

		// get all genes which are associated with the sub-trees
		GeneSet genes_children;
		foreach (const Phenotype& child, children)
		{
			genes_children.insert(db.phenotypeToGenes(db.phenotypeIdByAccession(child.accession()), true, false));
		}

		// intersect with genes present in the root node
		GeneSet genes_to_remove = genes_children.intersect(db.phenotypeToGenes(root_id, false, false));

		if (genes_to_remove.count() != 0)
		{
			// get phenotype id
			int pt_id = pt2id.value(root.accession());

			// remove all duplicate genes from parent node
			SqlQuery remove_gene_query = db.getQuery();
			remove_gene_query.prepare("DELETE FROM hpo_genes WHERE hpo_term_id=" + QByteArray::number(pt_id) + " AND gene=:0");
			foreach (const QByteArray& gene, genes_to_remove)
			{
				remove_gene_query.bindValue(0, gene);
				remove_gene_query.exec();
				removed_genes++;
			}
		}

		// start optimization for all child nodes
		foreach (const Phenotype& child, children)
		{
			optimizeHpoGeneTable(child, db, pt2id, removed_genes);
		}

		return;
	}
};

#include "main.moc"

int main(int argc, char *argv[])
{
	ConcreteTool tool(argc, argv);
	return tool.execute();
}
