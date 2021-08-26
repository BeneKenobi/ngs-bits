#include "ToolBase.h"
#include "Graph.h"
#include "OntologyTermCollection.h"
#include "Exceptions.h"
#include "Helper.h"
#include <QByteArray>
#include <QString>
#include <QHash>
#include <QList>
#include <QStringList>
#include <QTextStream>


struct NodeContent
{
    double no_of_diseases;
    double information_content;

    NodeContent()
        : no_of_diseases(0.0),
          information_content(0.0)
    {
    }
};

struct EdgeContent
{
    double weight;

    EdgeContent()
        : weight(0.0)
    {
    }
};

class ConcreteTool
        : public ToolBase
{
    Q_OBJECT

private:
    // create a directed graph from the ontology term collection
    Graph<NodeContent, EdgeContent> parseTermCollection(OntologyTermCollection &terms)
    {
        Graph<NodeContent, EdgeContent> hpo_graph(true);

        for(int i = 0; i < terms.size(); i++)
        {
            OntologyTerm term = terms.get(i);

            foreach(QByteArray parent, term.parentIDs())
            {
                NodeContent node_content_1;
                NodeContent node_content_2;
                EdgeContent edge_content;
                hpo_graph.addEdge(parent, node_content_1, term.id(), node_content_2, edge_content);
            }
        }

        return hpo_graph;
    }

    QHash<QString, QStringList> parseAnnotations(const QString &file)
    {
        QHash<QString, QStringList> annotations;

        QSharedPointer<QFile> reader = Helper::openFileForReading(file);
        QTextStream in(reader.data());

        while(!in.atEnd())
        {
            QString line = in.readLine();

            // skip comments
            if(line.startsWith("#"))
            {
                continue;
            }

            QStringList anno = line.split("\t");

            if(!annotations.contains(anno.at(3)))
            {
                annotations.insert(anno.at(3), QStringList());
            }

            if(!annotations[anno.at(3)].contains(anno.at(0)))
            {
                annotations[anno.at(3)].append(anno.at(0));
            }
        }

        return annotations;
    }

    QHash<QString, QStringList> parseAssociatedPhenotypes(const QString &file)
    {
        QHash<QString, QStringList> associated_terms;

        QSharedPointer<QFile> reader = Helper::openFileForReading(file);
        QTextStream in(reader.data());

        while(!in.atEnd())
        {
            QString line = in.readLine();

            // skip comments
            if(line.startsWith("#"))
            {
                continue;
            }

            QStringList entry = line.split("\t");

            if(!associated_terms.contains(entry.at(1)))
            {
                associated_terms.insert(entry.at(1), QStringList());
            }

            if(!associated_terms[entry.at(1)].contains(entry.at(2)))
            {
                associated_terms[entry.at(1)].append(entry.at(2));
            }
        }

        return associated_terms;
    }

public:
    ConcreteTool(int& argc, char *argv[])
        : ToolBase(argc, argv)
    {
    }

    virtual void setup()
    {
        setDescription("Creates simple representation of HPO directed acyclic graph.");
        addInfile("obo", "Human phenotype ontology obo file.", false);
        addInfile("annotation", "TSV file with HPO annotations to diseases", false);
        addInfile("associated-terms", "TSV file with all associated HPO terms for each gene", false);
        addOutfile("out", "Output TSV file with edges.", false);
    }

    virtual void main()
    {
        // init
        OntologyTermCollection terms(getInfile("obo"), true);
        Graph<NodeContent, EdgeContent> hpo_graph = parseTermCollection(terms);

        QHash<QString, QStringList> annotations = parseAnnotations(getInfile("annotation"));

        QHash<QString, QStringList> associated_terms = parseAssociatedPhenotypes(getInfile("associated-terms"));

        hpo_graph.store(getOutfile("out"));
    }
};

#include "main.moc"

int main(int argc, char *argv[])
{
    ConcreteTool tool(argc, argv);
    return tool.execute();
}

