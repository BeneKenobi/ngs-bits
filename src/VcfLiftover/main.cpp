#include "ToolBase.h"
#include "Settings.h"
#include "Exceptions.h"
#include "Helper.h"
#include "OntologyTermCollection.h"
#include "VcfFile.h"

#include <NGSHelper.h>

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
		setDescription("Lifts a VCF from GRCh37 to GRCh38.");

		addOutfile("failed", "Output VCF containig variants which cannot be lifted to GRCh38.", false, true);

		//optional
		addInfile("in", "Input VCF file. If unset, reads from STDIN.", true, true);
		addOutfile("out", "Output file. If unset, writes to STDOUT.", true, true);
		addInfile("ref_grch37", "GRCh37 reference genome FASTA file. If unset 'reference_genome' from the 'settings.ini' file is used.", true, false);
		addInfile("ref_grch38", "GRCh38 reference genome FASTA file. If unset 'reference_genome' from the 'settings.ini' file is used.", true, false);
		addFlag("debug", "Activates debug output.");

		changeLog(2021, 9, 21, "Initial implementation.");
	}

	virtual void main()
	{
		//init
		bool debug = getFlag("debug");
		QString in = getInfile("in");
		QString out = getOutfile("out");
		QString failed = getOutfile("failed");
		if(in!="" && in==out)
		{
			THROW(ArgumentException, "Input and output files must be different when streaming!");
		}

		// parse reference files
		QString ref_grch37 = getInfile("ref_grch37");
		if (ref_grch37 == "") ref_grch37 = Settings::string("reference_genome", true);
		if (ref_grch37 == "") THROW(CommandLineParsingException, "Reference genome FASTA for GRCh37 unset in both command-line and settings.ini file!");
		QString ref_grch38 = getInfile("ref_grch38");
		if (ref_grch38 == "") ref_grch38 = Settings::string("reference_genome_hg38", true);
		if (ref_grch38 == "") THROW(CommandLineParsingException, "Reference genome FASTA for GRCh38 unset in both command-line and settings.ini file!");



		//prepare output
		QSharedPointer<QFile> out_p = Helper::openFileForWriting(out, true);
		QTextStream out_stream(out_p.data());
		QSharedPointer<QFile> failed_p = Helper::openFileForWriting(failed, true);
		QTextStream failed_stream(failed_p.data());


		//load references
		QSharedPointer<FastaFileIndex> genome_index_grch37 = QSharedPointer<FastaFileIndex>(new FastaFileIndex(ref_grch37));
		QSharedPointer<FastaFileIndex> genome_index_grch38 = QSharedPointer<FastaFileIndex>(new FastaFileIndex(ref_grch38));

		//parse input
		QSharedPointer<QFile> in_p = Helper::openFileForReading(in, true);

		while(!in_p->atEnd())
		{
			QByteArray line = in_p->readLine();
			if (!line.endsWith('\n')) line += '\n';

			//skip empty lines
			if (line.trimmed().isEmpty()) continue;

			//write out headers/comments unchanged
			if (line.startsWith('#'))
			{
				out_p->write(line);
				failed_p->write(line);
				continue;
			}

			//parse variant
			QByteArrayList columns = line.split('\t');
			if (columns.count()<VcfFile::MIN_COLS)
			{
				THROW(FileParseException, "VCF line with too few columns in input file: \n" + line);
			}

			// parse position
			Chromosome chr = columns[VcfFile::CHROM];
			bool ok = false;
			int start = columns[VcfFile::POS].toInt(&ok);
			if (!ok)
			{
				THROW(FileParseException, "Could not convert VCF variant position '" + columns[1] + "' to integer!");
			}
			int end = start + columns[VcfFile::REF].length() - 1; //length of ref

			// parse sequences
			QByteArray ref = columns[VcfFile::REF];
			QByteArray obs = columns[VcfFile::ALT];

			//lift-over coordinates
			BedLine coords;
			try
			{
				coords = NGSHelper::liftOver(chr, start, end);
			}
			catch(Exception& e)
			{
				columns[VcfFile::INFO].append(";Error_Comment=Cannot_perform_liftover");
				failed_p->write(columns.join('\t'));
				if (debug) qDebug() << "Liftover failed for variant " << chr.str() << ":" << start << "-" << end << " " << ref << ">" << "obs ! (" << e.message() << ") ";
				continue;
			}

			//check new chromosome is ok
			if (!coords.chr().isNonSpecial())
			{
				columns[VcfFile::INFO].append(";Error_Comment=Liftover_on_special_chr");
				failed_p->write(columns.join('\t'));
				if (debug) qDebug() << "Liftover on special chr for variant " << chr.str() << ":" << start << "-" << end << " " << ref << ">" << "obs ! (" << coords.toString(true) << ") ";
				continue;
			}

			//check sequence context is the same (ref +-5 bases). If it is not, the strand might have changed, e.g. in NIPA1, GDF2, ANKRD35, TPTE, ...
			bool strand_changed = false;
			Sequence context_old = genome_index_grch37->seq(chr, start-5, 10 + ref.length());
			Sequence context_new = genome_index_grch38->seq(coords.chr(), coords.start()-5, 10 + ref.length());
			if (context_old!=context_new)
			{
				context_new.reverseComplement();
				if (context_old==context_new)
				{
					strand_changed = true;
				}
				else
				{
					context_new.reverseComplement();
					columns[VcfFile::INFO].append(";Error_Comment=Seq_context_changed(" + context_old + ">" + context_new + ")");
					failed_p->write(columns.join('\t'));
					if (debug) qDebug() << "Sequence context changed for variant " << chr.str() << ":" << start << "-" << end << " " << ref << ">" << "obs ! (" << context_old << ">" << context_new << ") ";
					continue;
				}
			}

			// valid liftover
			columns[VcfFile::CHROM] = coords.chr().strNormalized(true);
			columns[VcfFile::POS] = QByteArray::number(coords.start());
			out_p->write(columns.join('\t'));

		}


		out_p->flush();
		out_p->close();
		failed_p->flush();
		failed_p->close();
    }
};

#include "main.moc"

int main(int argc, char *argv[])
{
	ConcreteTool tool(argc, argv);
	return tool.execute();
}

