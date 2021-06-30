#ifndef PROCESSEDSAMPLEWIDGET_H
#define PROCESSEDSAMPLEWIDGET_H

#include "NGSD.h"

#include <QWidget>
#include <QLabel>
#include <QMenu>

namespace Ui {
class ProcessedSampleWidget;
}

class ProcessedSampleWidget
	: public QWidget
{
	Q_OBJECT


public:
	ProcessedSampleWidget(QWidget* parent, QString ps_id);
	~ProcessedSampleWidget();

	static void styleQualityLabel(QLabel* label, const QString& quality);

signals:
	void clearMainTableSomReport(QString ps_name);

protected slots:
	void updateGUI();
	void updateQCMetrics();
	void showPlot();
	void openSampleFolder();
	void openSampleTab();
	void openExternalDiseaseDatabase();
	void addRelation();
	void removeRelation();
	void editStudy();
	void addStudy();
	void removeStudy();
	void deleteSampleData();
	void loadVariantList();
	void queueSampleAnalysis();

	void openIgvTrack();
	void somRepDeleted();

	void openProcessedSampleTab(QString ps);
	void openRunTab(QString name);
	void openProjectTab(QString project_name);
	void openProcessingSystemTab(QString system_short_name);

	///Opens the processed sample edit dialog
	void edit();
	///Opens sample edit dialog
	void editSample();
	///Opens a dialog to edit the diagnostic status.
	void editDiagnosticStatus();
	///Opens a dialog to edit the disease group/info.
	void editDiseaseGroupAndInfo();
	///Opens a dialog to edit the disease details.
	void editDiseaseDetails();
	///Import sample relations for GenLab
	void importSampleRelations();

private:
	Ui::ProcessedSampleWidget* ui_;
	QString ps_id_;

	QString sampleName() const;
	QString processedSampleName() const;
	QString mergedSamples() const;
	void addIgvMenuEntry(QMenu* menu, PathType file_type);
};

#endif // PROCESSEDSAMPLEWIDGET_H
