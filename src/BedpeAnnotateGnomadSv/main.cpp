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

		changeLog(2020, 3, 20, "Initial commit.");
	}

	virtual void main()
	{
		//init
		QString in = getInfile("in");
		QString bed = getInfile("bed");
		QString out = getOutfile("out");



		//load annoation database
		BedFile gnomad_file;
		gnomad_file.load(bed);
		if (!gnomad_file.isSorted()) gnomad_file.sort();
		ChromosomalIndex<BedFile> anno_index(gnomad_file);



		//process BEDPE file
		BedpeFile bedpe_file;
		bedpe_file.load(in);

		// check if annotation already exisits:
		int i_gnomad = bedpe_file.annotationIndexByName("gnomAD", false);
		int i_gnomad_hom_hemi = bedpe_file.annotationIndexByName("gnomAD_hom_hemi", false);


		// create output buffer and copy comments and header
		QByteArrayList output_buffer;
		output_buffer.append(bedpe_file.comments());
		// get header
		QByteArrayList header = bedpe_file.annotationHeaders();
		// modify header if gene columns not already present
		if (i_gnomad < 0) header.append("gnomAD");
		if (i_gnomad_hom_hemi < 0) header.append("gnomAD_hom_hemi");
		// copy header
		output_buffer << "#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\t" + header.join("\t");


		for(int i=0; i<bedpe_file.count(); ++i)
		{
			BedpeLine line = bedpe_file[i];

			// find matching SV(s)			
			QVector<int> indices;
			if (line.type() == StructuralVariantType::BND)
			{
				indices = anno_index.matchingIndices(line.chr1(), line.start1(), line.end1());
			}
			else
			{
				indices = anno_index.matchingIndices(line.chr1(), line.start1(), line.end2());
			}

			foreach (int idx, indices)
			{
				// filter by type

				// exact match

				// extract gnomAD scores
			}


			QByteArray gnomad_af_string;
			QByteArray gnomad_hom_hemi_string;



//			BedFile affected_region = line.affectedRegion();

//			//determine annotations
//			QByteArrayList additional_annotations;

//			for(int j=0; j<affected_region.count(); ++j)
//			{
//				BedLine& region = affected_region[j];
//				QVector<int> indices = anno_index.matchingIndices(region.chr(), region.start(), region.end());
//				foreach(int index, indices)
//				{
//					const BedLine& match = gnomad_file[index];
//					bool anno_exists = match.annotations().count()>i_col;
//					if (anno_exists)
//					{
//						additional_annotations << match.annotations()[i_col];
//					}
//				}
//			}

			QList<QByteArray> annotations = line.annotations();
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
			line.setAnnotations(annotations);

			//add annotated line to buffer
			output_buffer << line.toTsv();
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
	}

};

#include "main.moc"

int main(int argc, char *argv[])
{
	ConcreteTool tool(argc, argv);
	return tool.execute();
}
