#include "ToolBase.h"
#include "Graph.h"
#include "OntologyTermCollection.h"
#include "Exceptions.h"
#include "Helper.h"
#include <QByteArray>
#include <QString>
#include <QHash>
#include <QSet>
#include <QList>
#include <QStringList>
#include <QTextStream>
#include <cmath>


struct NodeContent
{
    int associated_diseases;
    int directly_associated_diseases;
    double information_content;

    NodeContent()
        : associated_diseases(0),
          directly_associated_diseases(0),
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

    // parse diseases associated with each HPO term
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

            // file format: 1st column: DatabaseID; 4th column: HPO_ID
            QStringList anno = line.split("\t");

            if(!annotations.contains(anno.at(3)))
            {
                annotations.insert(anno.at(3), QStringList());
            }

            // add the disease (ID) to the list of associated diseases for the HPO term
            if(!annotations[anno.at(3)].contains(anno.at(0)))
            {
                annotations[anno.at(3)].append(anno.at(0));
            }
        }

        return annotations;
    }

    // parse HPO terms associated with each gene
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

            // file format: 2nd column: entrez-gene-symbol; 3rd column: HPO-Term-ID
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

    // write number of (directly) associated diseases into node content
    /*void countAssociatedDiseases(Graph<NodeContent, EdgeContent> &graph, QHash<QString, QStringList> &annotations)
    {
        Graph<NodeContent, EdgeContent>::NodePointer node;
        foreach(node, graph.adjacencyList().keys())
        {
            node.data()->nodeContent().directly_associated_diseases = annotations[node.data()->nodeName()].size();
        }
    }*/

    // recursively count number of directly and indirectly associated diseases for a term
    QSet<QString> countAssociatedDiseases(Graph<NodeContent, EdgeContent> &graph, QHash<QString, QStringList> &annotations,
                                QString node_name)
    {
        QSet<QString> associated_diseases = QSet<QString>::fromList(annotations[node_name]);
        graph.getNode(node_name).data()->nodeContent().directly_associated_diseases = associated_diseases.size();

        Graph<NodeContent, EdgeContent>::NodePointer node;
        foreach(node, graph.getAdjacentNodes(node_name))
        {
            associated_diseases.unite(countAssociatedDiseases(graph, annotations, node.data()->nodeName()));
        }

        graph.getNode(node_name).data()->nodeContent().associated_diseases = associated_diseases.size();

        return associated_diseases;
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

        // recursively count number of associated diseases, starting from root node "all"
        int total_diseases = countAssociatedDiseases(hpo_graph, annotations, "HP:0000001").size();

        // determine information content of each node
        Graph<NodeContent, EdgeContent>::NodePointer node;
        foreach(node, hpo_graph.adjacencyList().keys())
        {
            if(node.data()->nodeContent().associated_diseases != 0)
            {
                node.data()->nodeContent().information_content =
                        -log((double) node.data()->nodeContent().associated_diseases / total_diseases);
            }
            else
            {
                node.data()->nodeContent().information_content = 0;
            }
        }

        /*QTextStream out(stdout);
        //Graph<NodeContent, EdgeContent>::NodePointer node;
        foreach(node, hpo_graph.adjacencyList().keys())
        {
            out << node.data()->nodeName() << "\t"
                << node.data()->nodeContent().directly_associated_diseases << "\t"
                << node.data()->nodeContent().associated_diseases << "\t"
                << node.data()->nodeContent().information_content << endl;
        }
        out << total_diseases << endl;*/

        hpo_graph.store(getOutfile("out"));
    }
};

#include "main.moc"

int main(int argc, char *argv[])
{
    ConcreteTool tool(argc, argv);
    return tool.execute();
}

