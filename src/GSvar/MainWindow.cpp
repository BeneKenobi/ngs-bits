#include "MainWindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include "Settings.h"
#include "Exceptions.h"
#include "ChromosomalIndex.h"
#include "Log.h"
#include "Helper.h"
#include "GUIHelper.h"
#include "GeneSet.h"
#include <QDir>
#include <QBitArray>
#include <QDesktopServices>
#include <QUrl>
#include <QTcpSocket>
#include <QTime>
#include <ProxyDataService.h>
#include "ExternalToolDialog.h"
#include "ReportDialog.h"
#include <QBrush>
#include <QFont>
#include <QInputDialog>
#include <QClipboard>
#include <QProgressBar>
#include <QToolButton>
#include <QMimeData>
#include <QSqlError>
#include <QChartView>
#include <GenLabDB.h>
#include <QToolTip>
#include <QProcess>
#include <QImage>
#include <QBuffer>
QT_CHARTS_USE_NAMESPACE
#include "ReportWorker.h"
#include "ScrollableTextDialog.h"
#include "AnalysisStatusWidget.h"
#include "HttpHandler.h"
#include "ValidationDialog.h"
#include "ClassificationDialog.h"
#include "BasicStatistics.h"
#include "ApprovedGenesDialog.h"
#include "GeneWidget.h"
#include "PhenoToGenesDialog.h"
#include "GenesToRegionsDialog.h"
#include "SubpanelDesignDialog.h"
#include "SubpanelArchiveDialog.h"
#include "GapDialog.h"
#include "EmailDialog.h"
#include "CnvWidget.h"
#include "CnvList.h"
#include "RohWidget.h"
#include "GeneSelectorDialog.h"
#include "NGSHelper.h"
#include "QCCollection.h"
#include "DiseaseInfoWidget.h"
#include "SmallVariantSearchWidget.h"
#include "TSVFileStream.h"
#include "OntologyTermCollection.h"
#include "SvWidget.h"
#include "VariantWidget.h"
#include "SomaticReportConfigurationWidget.h"
#include "SingleSampleAnalysisDialog.h"
#include "MultiSampleDialog.h"
#include "TrioDialog.h"
#include "SomaticDialog.h"
#include "Histogram.h"
#include "ProcessedSampleWidget.h"
#include "DBSelector.h"
#include "SequencingRunWidget.h"
#include "SimpleCrypt.h"
#include "ToolBase.h"
#include "BedpeFile.h"
#include "SampleSearchWidget.h"
#include "ProcessedSampleSelector.h"
#include "ReportVariantDialog.h"
#include "SomaticReportVariantDialog.h"
#include "GSvarHelper.h"
#include "SampleDiseaseInfoWidget.h"
#include "QrCodeFactory.h"
#include "SomaticRnaReport.h"
#include "ProcessingSystemWidget.h"
#include "ProjectWidget.h"
#include "DBEditor.h"
#include "TsvTableWidget.h"
#include "DBTableAdministration.h"
#include "SequencingRunOverview.h"
#include "MidCheckWidget.h"
#include "CnvSearchWidget.h"
#include "VariantValidationWidget.h"
#include "SomaticReportDialog.h"
#include "GeneOmimInfoWidget.h"
#include "LoginManager.h"
#include "LoginDialog.h"
#include "IgvSessionManager.h"
#include "GeneInfoDBs.h"
#include "VariantConversionWidget.h"
#include "PasswordDialog.h"
#include "CircosPlotWidget.h"
#include "SomaticReportSettings.h"
#include "CytobandToRegionsDialog.h"
#include "RepeatExpansionWidget.h"
#include "SomaticDataTransferWidget.h"
#include "PRSWidget.h"
#include "EvaluationSheetEditDialog.h"
#include "SvSearchWidget.h"
#include "PublishedVariantsWidget.h"
#include "PreferredTranscriptsWidget.h"
#include "TumorOnlyReportWorker.h"
#include "TumorOnlyReportDialog.h"
#include "VariantScores.h"
#include "CfDNAPanelDesignDialog.h"
#include "DiseaseCourseWidget.h"
#include "CfDNAPanelWidget.h"
#include "SomaticVariantInterpreterWidget.h"
#include "AlleleBalanceCalculator.h"
#include "ExpressionGeneWidget.h"
#include "GapClosingDialog.h"
#include "XmlHelper.h"
#include "GermlineReportGenerator.h"
#include "SomaticReportHelper.h"
#include "Statistics.h"
#include "NGSDReplicationWidget.h"
#include "CohortAnalysisWidget.h"
#include "cfDNARemovedRegions.h"
#include "CfDNAPanelBatchImport.h"
#include "ClinvarUploadDialog.h"
#include "GenomeVisualizationWidget.h"
#include "LiftOverWidget.h"
#include "CacheInitWorker.h"
#include "BlatWidget.h"
#include "FusionWidget.h"
#include "CohortExpressionDataWidget.h"
#include "CausalVariantEditDialog.h"
#include "VariantOpenDialog.h"
#include "GeneSelectionDialog.h"
#include "ExpressionOverviewWidget.h"
#include "ExpressionExonWidget.h"
#include "SplicingWidget.h"
#include "VariantHgvsAnnotator.h"
#include "VirusDetectionWidget.h"
#include "SomaticcfDNAReport.h"
#include "MaintenanceDialog.h"
#include "ClientHelper.h"
#include "ProxyDataService.h"
#include "RefGenomeService.h"
#include "GHGAUploadDialog.h"
#include "BurdenTestWidget.h"
#include "IgvLogWidget.h"
#include "SettingsDialog.h"
#include "GlobalServiceProvider.h"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui_()
	, var_last_(-1)
	, busy_dialog_(nullptr)
	, notification_label_(new QLabel())
    , igv_history_label_(new ClickableLabel())
	, filename_()
	, variants_changed_()
	, last_report_path_(QDir::homePath())
    , init_timer_(this, true)
    , server_version_()
{
	//setup GUI
	ui_.setupUi(this);
	setWindowTitle(appName());
	GUIHelper::styleSplitter(ui_.splitter);
	ui_.splitter->setStretchFactor(0, 10);
	ui_.splitter->setStretchFactor(1, 1);
	GUIHelper::styleSplitter(ui_.splitter_2);
	ui_.splitter_2->setStretchFactor(0, 10);
	ui_.splitter_2->setStretchFactor(1, 1);
	connect(ui_.tabs, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTab(int)));
	ui_.actionDebug->setVisible(Settings::boolean("debug_mode_enabled", true));
	ui_.actionEncrypt->setVisible(Settings::boolean("debug_mode_enabled", true));

	// add rna menu
	rna_menu_btn_ = new QToolButton();
	rna_menu_btn_->setObjectName("rna_btn");
	rna_menu_btn_->setIcon(QIcon(":/Icons/RNA.png"));
	rna_menu_btn_->setToolTip("Open RNA menu entries");
	rna_menu_btn_->setMenu(new QMenu());
	rna_menu_btn_->menu()->addAction(ui_.actionExpressionData);
	rna_menu_btn_->menu()->addAction(ui_.actionExonExpressionData);
	rna_menu_btn_->menu()->addAction(ui_.actionShowSplicing);
	rna_menu_btn_->menu()->addAction(ui_.actionShowRnaFusions);
	rna_menu_btn_->menu()->addAction(ui_.actionShowProcessingSystemCoverage);
	rna_menu_btn_->setPopupMode(QToolButton::InstantPopup);

	ui_.actionExpressionData->setEnabled(false);
	ui_.actionExonExpressionData->setEnabled(false);
	ui_.actionShowSplicing->setEnabled(false);
	ui_.actionShowRnaFusions->setEnabled(false);

	ui_.tools->addWidget(rna_menu_btn_);


	// add cfdna menu
	cfdna_menu_btn_ = new QToolButton();
	cfdna_menu_btn_->setObjectName("cfdna_btn");
	cfdna_menu_btn_->setIcon(QIcon(":/Icons/cfDNA.png"));
	cfdna_menu_btn_->setToolTip("Open cfDNA menu entries");
	cfdna_menu_btn_->setMenu(new QMenu());
	cfdna_menu_btn_->menu()->addAction(ui_.actionDesignCfDNAPanel);
	cfdna_menu_btn_->menu()->addAction(ui_.actionShowCfDNAPanel);
	cfdna_menu_btn_->menu()->addAction(ui_.actionCfDNADiseaseCourse);
	cfdna_menu_btn_->menu()->addAction(ui_.actionCfDNAAddExcludedRegions);
	cfdna_menu_btn_->setPopupMode(QToolButton::InstantPopup);
	ui_.tools->addWidget(cfdna_menu_btn_);
	// deaktivate on default (only available in somatic)
	cfdna_menu_btn_->setVisible(false);
	cfdna_menu_btn_->setEnabled(false);

	//signals and slots
    connect(ui_.actionExit, SIGNAL(triggered()), this, SLOT(closeAndLogout()));

	connect(ui_.filters, SIGNAL(filtersChanged()), this, SLOT(refreshVariantTable()));
	connect(ui_.vars, SIGNAL(itemSelectionChanged()), this, SLOT(updateVariantDetails()));
	connect(ui_.vars, SIGNAL(cellDoubleClicked(int, int)), this, SLOT(variantCellDoubleClicked(int, int)));
	connect(ui_.vars, SIGNAL(showMatchingCnvsAndSvs(BedLine)), this, SLOT(showMatchingCnvsAndSvs(BedLine)));
	connect(ui_.vars->verticalHeader(), SIGNAL(sectionDoubleClicked(int)), this, SLOT(variantHeaderDoubleClicked(int)));
	ui_.vars->verticalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(ui_.vars->verticalHeader(), SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(varHeaderContextMenu(QPoint)));

	connect(ui_.actionDesignSubpanel, SIGNAL(triggered()), this, SLOT(openSubpanelDesignDialog()));
	connect(ui_.filters, SIGNAL(phenotypeImportNGSDRequested()), this, SLOT(importPhenotypesFromNGSD()));
	connect(ui_.filters, SIGNAL(phenotypeSubPanelRequested()), this, SLOT(createSubPanelFromPhenotypeFilter()));

	//variants tool bar
	connect(ui_.vars_copy_btn, SIGNAL(clicked(bool)), ui_.vars, SLOT(copyToClipboard()));
	connect(ui_.vars_resize_btn, SIGNAL(clicked(bool)), ui_.vars, SLOT(adaptColumnWidthsCustom()));
	ui_.vars_export_btn->setMenu(new QMenu());
	ui_.vars_export_btn->menu()->addAction("Export GSvar (filtered)", this, SLOT(exportGSvar()));
	ui_.vars_export_btn->menu()->addAction("Export VCF (filtered)", this, SLOT(exportVCF()));
	ui_.report_btn->setMenu(new QMenu());
	ui_.report_btn->menu()->addAction(QIcon(":/Icons/Report_add_causal.png"), "Add/edit other causal variant", this, SLOT(editOtherCausalVariant()));
	ui_.report_btn->menu()->addAction(QIcon(":/Icons/Report_exclude.png"), "Delete other causal variant", this, SLOT(deleteOtherCausalVariant()));
	ui_.report_btn->menu()->addSeparator();
	ui_.report_btn->menu()->addAction(QIcon(":/Icons/Report.png"), "Generate report", this, SLOT(generateReport()));
	ui_.report_btn->menu()->addAction(QIcon(":/Icons/Report.png"), "Generate evaluation sheet", this, SLOT(generateEvaluationSheet()));
	ui_.report_btn->menu()->addAction(QIcon(":/Icons/Report_info.png"), "Show report configuration info", this, SLOT(showReportConfigInfo()));
	ui_.report_btn->menu()->addSeparator();
	ui_.report_btn->menu()->addAction(QIcon(":/Icons/Report_finalize.png"), "Finalize report configuration", this, SLOT(finalizeReportConfig()));
	ui_.report_btn->menu()->addSeparator();
	ui_.report_btn->menu()->addAction("Transfer somatic data to MTB", this, SLOT(transferSomaticData()) );
	connect(ui_.vars_folder_btn, SIGNAL(clicked(bool)), this, SLOT(openVariantListFolder()));
	connect(ui_.open_qc_files, SIGNAL(clicked(bool)), this, SLOT(openVariantListQcFiles()));
	ui_.vars_ranking->setMenu(new QMenu());
	ui_.vars_ranking->menu()->addAction("dominant model", this, SLOT(variantRanking()))->setObjectName("GSvar_v2_dominant");
	ui_.vars_ranking->menu()->addAction("recessive model", this, SLOT(variantRanking()))->setObjectName("GSvar_v2_recessive");
	ui_.vars_af_hist->setMenu(new QMenu());
	ui_.vars_af_hist->menu()->addAction("Show AF histogram (all small variants)", this, SLOT(showAfHistogram_all()));
	ui_.vars_af_hist->menu()->addAction("Show AF histogram (small variants after filter)", this, SLOT(showAfHistogram_filtered()));
	ui_.vars_af_hist->menu()->addSeparator();
	ui_.vars_af_hist->menu()->addAction("Show CN histogram (in given region)", this, SLOT(showCnHistogram()));
	ui_.vars_af_hist->menu()->addAction("Show BAF histogram (in given region)", this, SLOT(showBafHistogram()));

	connect(ui_.ps_details, SIGNAL(clicked(bool)), this, SLOT(openProcessedSampleTabsCurrentAnalysis()));

	//if at home, use Patientenserver
	QString gsvar_report_folder = Settings::path("gsvar_report_folder", true);
	if (gsvar_report_folder!="" && QDir(gsvar_report_folder).exists())
	{
		last_report_path_ = gsvar_report_folder;
	}

	//add notification icon
	notification_label_->hide();
	notification_label_->setScaledContents(true);
	notification_label_->setMaximumSize(16,16);
	notification_label_->setPixmap(QPixmap(":/Icons/email.png"));
	ui_.statusBar->addPermanentWidget(notification_label_);	
	ui_.actionVirusDetection->setEnabled(false);

    igv_history_label_->setScaledContents(true);
    igv_history_label_->setMaximumSize(16,16);
    igv_history_label_->setPixmap(QPixmap(":/Icons/IGV.png"));
    igv_history_label_->setToolTip("Show the history of IGV commands");
    ui_.statusBar->addPermanentWidget(igv_history_label_);
    connect(igv_history_label_, SIGNAL(clicked(QPoint)), this, SLOT(displayIgvHistoryTable(QPoint)));

	//init cache in background thread (it takes about 6 seconds)
	CacheInitWorker* worker = new CacheInitWorker();
	worker->start();

	// Setting a value for the current working directory. On Linux it is defined in the TMPDIR environment
	// variable or /tmp if TMPDIR is not set. On Windows it is saved in the TEMP or TMP environment variable.
	// e.g. c:\Users\USER_NAME\AppData\Local\Temp
	// It is needed to enable saving *.bai files while accessing remote *.bam files. htsLib tries to save the
	// index file locally, if it deals with a remote *.bam file. The index file is always saved at the current
	// working directory, and it seems there is no way to change it. On some systems users may not have write
	// priveleges for the working directory and this is precisely why we came up with this workaround:
	QDir::setCurrent(QDir::tempPath());

	//enable timers needed in client-server mode
	if (ClientHelper::isClientServerMode())
	{
		// renew existing session, if it is about to expire
		// a new token will be requested slightly in advance
		QTimer *login_timer = new QTimer(this);
		connect(login_timer, &QTimer::timeout, this, &LoginManager::renewLogin);
		login_timer->start(20 * 60 * 1000); // every 20 minutes

		//check if the server is running
		QTimer *server_ping_timer = new QTimer(this);
		connect(server_ping_timer, SIGNAL(timeout()), this, SLOT(checkServerAvailability()));
		server_ping_timer->start(10 * 60 * 1000); // every 10 minutes


		//check if there are new notifications for the users
		if (Settings::boolean("display_user_notifications", true))
		{
			QTimer *user_notification_timer = new QTimer(this);
			connect(user_notification_timer, SIGNAL(timeout()), this, SLOT(checkUserNotifications()));
			user_notification_timer->start(12 * 60 * 1000); // every 12 minutes
		}

		displayed_maintenance_message_id_ = "";
	}

	connect(ui_.vars, SIGNAL(publishToClinvarTriggered(int, int)), this, SLOT(uploadToClinvar(int, int)));
	connect(ui_.vars, SIGNAL(alamutTriggered(QAction*)), this, SLOT(openAlamut(QAction*)));

	// Environment variable containing the file path to the list of certificate authorities
	// (needed for HTTPS to work correctly, especially for htslib and BamReader)
	QString curl_ca_bundle = Settings::string("curl_ca_bundle", true);
	if ((Helper::isWindows()) && (!curl_ca_bundle.isEmpty()))
	{
		if (!qputenv("CURL_CA_BUNDLE", curl_ca_bundle.toUtf8()))
		{
			Log::error("Could not set CURL_CA_BUNDLE variable, access to BAM files over HTTPS may not be possible");
		}
	}

	update_info_toolbar_ = new QToolBar;
	update_info_toolbar_->hide();
	addToolBar(Qt::TopToolBarArea, update_info_toolbar_);
}

QString MainWindow::appName() const
{
	QString name = QCoreApplication::applicationName();

	GenomeBuild build = GSvarHelper::build();
	if (build!=GenomeBuild::HG38) name += " - " + buildToString(build);

	return name;
}

bool MainWindow::isServerRunning()
{
    int status_code = -1;
    ServerInfo server_info = ClientHelper::getServerInfo(status_code);

    if (server_info.isEmpty())
	{
		QMessageBox::warning(this, "Server not available", "GSvar is configured for the client-server mode, but the server is not available. The application will be closed");
		return false;
	}

    if (status_code!=200)
    {
        QMessageBox::warning(this, "Server availability problem", "Server replied with " + QString::number(status_code) + " code. The application will be closed");
        return false;
    }

    if (!server_version_.isEmpty() && (server_version_ != server_info.version))
    {
        QMessageBox::information(this, "Server version changed", "Server version has changed from " + server_version_ + " to " + server_info.version + ". No action is required");
    }
    server_version_ = server_info.version;

	if (ClientHelper::serverApiVersion() != server_info.api_version)
	{
		QMessageBox::warning(this, "Version mismatch", "GSvar uses API " + ClientHelper::serverApiVersion() + ", while the server uses API " + server_info.api_version + ". No stable work can be guaranteed. The application will be closed");
		return false;
	}

    return true;
}

void MainWindow::checkServerAvailability()
{
	if (!isServerRunning())
	{
		close();
	}
}

void MainWindow::checkUserNotifications()
{
	UserNotification user_notification = ClientHelper::getUserNotification();

	if (user_notification.id.isEmpty() || user_notification.message.isEmpty()) return;
	if ((!displayed_maintenance_message_id_.isEmpty()) && (displayed_maintenance_message_id_ == user_notification.id)) return;

	displayed_maintenance_message_id_ = user_notification.id;
	QMessageBox::information(this, "Important information", user_notification.message);
}

void MainWindow::checkClientUpdates()
{
	if (!ClientHelper::isClientServerMode()) return;

	ClientInfo client_info = ClientHelper::getClientInfo();
	if (client_info.isEmpty())
	{
		Log::warn("Could not retrieve updates information from the server");
		return;
	}

	int commit_pos = QCoreApplication::applicationVersion().lastIndexOf("-");
	if (commit_pos>-1)
	{
		QString short_version = QCoreApplication::applicationVersion().left(commit_pos);
		if (ClientInfo(short_version, "").isOlderThan(client_info))
		{
			Log::info("Client version from the server: " + client_info.version);
			update_info_toolbar_->clear();
			update_info_toolbar_->setAutoFillBackground(true);
			update_info_toolbar_->setStyleSheet("QToolBar {background: red;}");
			QLabel *update_info_label = new QLabel;
			update_info_label->setText(client_info.message.isEmpty() ? "Please restart the application" : client_info.message);
			update_info_label->setStyleSheet("QLabel {color: white}");
			update_info_toolbar_->addWidget(update_info_label);
			update_info_toolbar_->show();
		}
		else
		{
			update_info_toolbar_->hide();
        }
    }
}

QString MainWindow::getCurrentFileName()
{
    return filename_;
}

AnalysisType MainWindow::getCurrentAnalysisType()
{
    return variants_.type();
}

void MainWindow::on_actionDebug_triggered()
{
	QTime timer;
	timer.start();

	QString user = Helper::userName();
	qDebug() << user;
	if (user=="ahsturm1")
	{
		//VariantHgvsAnnotator debugging
		/*
		QString genome_file = Settings::string("reference_genome", false);
		FastaFileIndex genome_idx(genome_file);
		VariantHgvsAnnotator hgvs_anno(genome_idx);

		VcfLine variant = VcfLine("chr7", 157009948, "CA", QList<Sequence>() << "CCGCGGCGGCG");
		qDebug() << variant.toString(true);
		TranscriptList transcripts = NGSD().transcriptsOverlapping(variant.chr(), variant.start(), variant.end());
		foreach(const Transcript& t, transcripts)
		{
			qDebug() << t.name();
			VariantConsequence hgvs = hgvs_anno.annotate(t, variant, true);
			qDebug() << hgvs.toString();
		}
		*/

		//extract VCF with variants that have class 4/5
		/*
		VariantList variants;
		NGSD db;
		QList<int> variant_ids = db.getValuesInt("SELECT DISTINCT v.id FROM `variant` v, variant_classification vc WHERE vc.variant_id= v.id AND (vc.class='4' OR vc.class='5')");
		foreach(int var_id, variant_ids)
		{
			Variant v = db.variant(QString::number(var_id));
			variants.append(v);
		}

		VcfFile vcf = VcfFile::fromGSvar(variants, genome_file);
		vcf.sort(false);
		vcf.store("C:\\Marc\\class4_and_5.vcf");
		*/

		//export a transcript definition from NGSD for VariantHgvsAnnotator test
		/*
		QByteArray trans = "ENST00000252971";
		NGSD db;
		int trans_id = db.transcriptId(trans);
		Transcript t = db.transcript(trans_id);

		QTextStream out(stdout);
		out << "Transcript trans_" << t.gene() << "()" << endl;
		out << "{" << endl;
		out << "\tTranscript t;" << endl;
		out << "\tt.setGene(\"" << t.gene() << "\");" << endl;
		out << "\tt.setName(\"" << t.name() << "\");" << endl;
		out << "\tt.setVersion(" << t.version() << ");" << endl;
		out << "\tt.setSource(Transcript::ENSEMBL);" << endl;
		out << "\tt.setStrand(Transcript::" << (t.isPlusStrand() ? "PLUS" : "MINUS") << ");" << endl;
		out << "\t" << endl;
		out << "\tBedFile regions;" << endl;
		for (int i=0; i<t.regions().count(); ++i)
		{
			const BedLine& line = t.regions()[i];
			out << "\tregions.append(BedLine(\"" << line.chr().str() << "\", " << line.start() << ", " << line.end() << "));" << endl;
		}
		out << "\tt.setRegions(regions";
		if (t.isCoding())
		{
			out << ", " << t.codingStart();
			out << ", " << t.codingEnd();
		}
		out << ");" << endl;
		out << "\t" << endl;
		out << "\treturn t;" << endl;
		out << "}" << endl;
		*/

		//generate somatic XML files in batch
		/*
		NGSD db;
		QStringList pairs;
		pairs << "DNA2203823A1_01	DNA2202978A1_01";
		foreach(QString pair, pairs)
		{
			QStringList parts =  pair.split("\t");
			QString ps_tumor = parts[0];
			QString ps_normal = parts[1];
			QStringList errors;
			QString xml_file = "C:\\Marc\\somatic_reports\\"+ps_tumor+"-"+ps_normal+".xml";

			try
			{
				QString ps_tumor_id = db.processedSampleId(ps_tumor);
				QString ps_normal_id = db.processedSampleId(ps_normal);

				QStringList somatic_analyses =  GlobalServiceProvider::database().secondaryAnalyses(ps_tumor + "-" + ps_normal, "somatic");
				if (somatic_analyses.isEmpty())
				{
					THROW(Exception, "Could not find secondary analysis in NGSD");
				}

				QString gsvar = somatic_analyses[0];
				loadFile(gsvar, true);
				ui_.filters->setFilter("somatic");
				ui_.filters->setTargetRegionByDisplayName("Sure Select Somatic Cancer Panel v5");
				somatic_report_settings_.target_region_filter = ui_.filters->targetRegion();

				somatic_report_settings_.filters = ui_.filters->filters();
				somatic_report_settings_.preferred_transcripts = GSvarHelper::preferredTranscripts();

				//get SO-terms
				OntologyTermCollection obo_terms("://Resources/so-xp_3_1_0.obo", true);
				QList<QByteArray> ids;
				ids << obo_terms.childIDs("SO:0001580",true); //coding variants
				ids << obo_terms.childIDs("SO:0001568",true); //splicing variants
				foreach(const QByteArray& id, ids)
				{
					somatic_report_settings_.obo_terms_coding_splicing.add(obo_terms.getByID(id));
				}

				//get phenotype infos from NGSD
				QStringList tmp_icd10;
				QStringList tmp_phenotype;
				foreach(const auto& entry, db.getSampleDiseaseInfo(db.sampleId(ps_tumor)) )
				{
					if(entry.type == "ICD10 code") tmp_icd10.append(entry.disease_info);
					if(entry.type == "clinical phenotype (free text)") tmp_phenotype.append(entry.disease_info);
				}
				somatic_report_settings_.icd10 = tmp_icd10.join(", ");
				somatic_report_settings_.phenotype = tmp_phenotype.join(", ");

				//load IGV screenshot
				if(GlobalServiceProvider::fileLocationProvider().getSomaticIgvScreenshotFile().exists)
				{
					QImage picture;
					picture = QImage(GlobalServiceProvider::fileLocationProvider().getSomaticIgvScreenshotFile().filename);
					if(!picture.isNull())
					{
						if( (uint)picture.width() > 1200 ) picture = picture.scaledToWidth(1200, Qt::TransformationMode::SmoothTransformation);
						if( (uint)picture.height() > 1200 ) picture = picture.scaledToHeight(1200, Qt::TransformationMode::SmoothTransformation);

						QByteArray png_data = "";
						QBuffer buffer(&png_data);
						buffer.open(QIODevice::WriteOnly);
						if(picture.save(&buffer, "PNG"))
						{
							somatic_report_settings_.igv_snapshot_png_hex_image = png_data.toHex();
							somatic_report_settings_.igv_snapshot_width = picture.width();
							somatic_report_settings_.igv_snapshot_height = picture.height();
						}
					}
				}

				//load germline variants
				VariantList variants_germline;
				variants_germline.load(GlobalServiceProvider::database().processedSamplePath(ps_normal_id, PathType::GSVAR).filename);
				if(!SomaticReportHelper::checkGermlineSNVFile(variants_germline))
				{
					THROW(Exception, "DNA report cannot be created because germline GSVar file is invalid. Please check control tissue variant file.");
				}

				QStringList messages;
				somatic_report_settings_.report_config = db.somaticReportConfig(ps_tumor_id, ps_normal_id, variants_, cnvs_, somatic_control_tissue_variants_, messages);
				if(!messages.isEmpty())
				{
					THROW(Exception, "Report config problems: " + messages.join("\n"));
				}

				//Store XML file with the same somatic report configuration settings
				SomaticReportHelper report(GSvarHelper::build(), variants_, cnvs_, variants_germline, somatic_report_settings_);
				report.storeXML(xml_file);
			}
			catch(Exception e)
			{
				errors << e.message();
			}

			//log
			QString log_file = "C:\\Marc\\somatic_reports\\logs\\"+ps_tumor+"-"+ps_normal+".log";
			if (QFile::exists(log_file)) QFile::remove(log_file);
			if (errors.count()>0)
			{
				Helper::storeTextFile(log_file, errors);
			}
		}
		*/

		//show genes without transcripts
		/*
		NGSD db;
		QSet<int> gene_ids = db.getValuesInt("SELECT id FROM gene").toSet();
		qDebug() << gene_ids.count();
		QSet<int> gene_ids_with_transcript = db.getValuesInt("SELECT DISTINCT gene_id FROM gene_transcript").toSet();
		qDebug() << gene_ids_with_transcript.count();
		gene_ids.subtract(gene_ids_with_transcript);
		qDebug() << gene_ids.count();
		foreach(int gene_id, gene_ids)
		{
			QString type = db.getValue("SELECT type FROM gene WHERE id="+QString::number(gene_id)).toString();
			if (type=="protein-coding gene") QTextStream(stdout) << db.getValue("SELECT symbol FROM gene WHERE id="+QString::number(gene_id)).toString() << "\t" << type << endl;
		}
		*/

		//Delete report config CNVs of samples that where the report configuration was not changed since 06.12.22 (for re-import of CNV report config data from HG19 databases - necessary because of CNV calling bug at chromosome ends)
		/*
		NGSD db;
		QList<int> rc_ids_with_cnv_rc = db.getValuesInt("SELECT DISTINCT rc.id FROM report_configuration rc, report_configuration_cnv rcc WHERE rc.id=rcc.report_configuration_id AND rc.last_edit_date < \"2021-12-06\" AND rc.created_date < \"2021-12-06\"");
		foreach(int rc_id, rc_ids_with_cnv_rc)
		{
			QString ps_id = db.getValue("SELECT processed_sample_id FROM report_configuration WHERE id=:0", false, QString::number(rc_id)).toString();
			qDebug() << "Deleting report config CNVs of " << db.processedSampleName(ps_id) << "ps_id=" << ps_id  << "rc_id=" << rc_id;
			SqlQuery query = db.getQuery();
			query.exec("DELETE FROM `report_configuration_cnv` WHERE `report_configuration_id`='"+QString::number(rc_id)+"'");
			qDebug() << "  Affected rows:" << query.numRowsAffected();
		}
		*/

		//Delete small variants report config of samples that have no variants imported (caused by error in NGSDReplicationWidget)
		/*
		NGSD db;
		QList<int> ps_ids_with_small_variant_rc = db.getValuesInt("SELECT DISTINCT rc.processed_sample_id FROM report_configuration rc, report_configuration_variant rcv WHERE rc.id=rcv.report_configuration_id");
		qDebug() << ps_ids_with_small_variant_rc.count();
		QSet<int> ps_ids_with_variants_imported = db.getValuesInt("SELECT DISTINCT(processed_sample_id) FROM `detected_variant`").toSet();
		qDebug() << ps_ids_with_variants_imported.count();
		foreach(int ps_id, ps_ids_with_small_variant_rc)
		{
			if (!ps_ids_with_variants_imported.contains(ps_id))
			{
				QString ps_id_str = QString::number(ps_id);
				QString rc_id = db.reportConfigId(ps_id_str);
				qDebug() << "Deleting " << db.processedSampleName(ps_id_str) << "ps_id=" << ps_id_str  << "rc_id=" << rc_id;
				db.getQuery().exec("DELETE FROM `report_configuration_variant` WHERE `report_configuration_id`='"+rc_id+"'");
			}
		}
		*/

		//Check HPO terms in NGSD
		/*
		NGSD db;
		PhenotypeList valid_terms;
		valid_terms << db.phenotypeChildTerms(db.phenotypeIdByName("Clinical course"), true);
		valid_terms << db.phenotypeChildTerms(db.phenotypeIdByName("Clinical modifier"), true);
		valid_terms << db.phenotypeChildTerms(db.phenotypeIdByName("Past medical history"), true);
		valid_terms << db.phenotypeChildTerms(db.phenotypeIdByName("Phenotypic abnormality"), true);

		auto file = Helper::openFileForWriting("C:\\Marc\\hpos.tsv");
		QTextStream stream(file.data());
		stream << "#hpo_id\thpo_name\tsamples\terrors\n";
		SqlQuery query = db.getQuery();
		query.exec("SELECT `disease_info`, COUNT(`sample_id`) as sample_count FROM `sample_disease_info` WHERE `type`='HPO term id' GROUP BY disease_info");
		while(query.next())
		{
			QStringList errors;

			QByteArray hpo_id = query.value(0).toByteArray().trimmed();
			QString samples = query.value(1).toString().trimmed();

			QString hpo_name;
			try
			{
				int id = db.phenotypeIdByAccession(hpo_id);

				hpo_name = db.phenotype(id).name();

				if (!valid_terms.containsAccession(hpo_id))
				{
					errors << "Not a child of 'phenotypic abnormality', 'Clinical modifier' or 'Past medical history'";
				}
			}
			catch(Exception & e)
			{
				errors << e.message();
			}

			stream << hpo_id << "\t"  << hpo_name << "\t"  << samples << "\t" << errors.join(", ") << "\n";
			stream.flush();
		}
		*/

		//export of recurring variants with similar phenotype
		/*
		NGSD db;
		auto file = Helper::openFileForWriting("C:\\Marc\\vars_"+Helper::dateTime("").replace(":", "")+".tsv");
		QTextStream stream(file.data());
		stream << "#gene\ttranscript\tvariant\tHGVS.p\ttype\timpact\tgnomad_AF\tclassification\tnum_affected\tnum_unaffeacted\tnum_unknown\tshared_disease_group\tsamples_with_hpo\tshared_hpo_term\n";

		//NGSD db;
		int c_gene = 0;
		QStringList genes = db.getValues("SELECT DISTINCT gene FROM omim_gene og WHERE id IN (SELECT DISTINCT omim_gene_id  FROM omim_phenotype) ORDER BY gene ASC");
		foreach(QString gene, genes)
		{
			qDebug() << ++c_gene << gene;
			int gene_id = db.geneToApprovedID(gene.toUtf8());
			if (gene_id==-1)
			{
				stream << "##" << gene << ": skipped - no approved gene name\n";
				continue;
			}

			Transcript lct = db.longestCodingTranscript(gene_id, Transcript::ENSEMBL, false, false);
			if (lct.codingRegions().baseCount()==0)
			{
				stream << "##" << gene << ": skipped - no longest coding transcript\n";
				continue;
			}

			BedFile roi_coding = lct.codingRegions();
			roi_coding.extend(20);
			roi_coding.merge();


			SqlQuery query = db.getQuery();
			QString af = "0.001";
			Chromosome chr = roi_coding[0].chr();
			query.exec("SELECT v.id, v.start, v.end, v.ref, v.obs, v.coding, v.gnomad FROM variant v WHERE chr='" + chr.strNormalized(true)  + "' AND start>='" + QString::number(roi_coding[0].start()) + "' AND end<='" + QString::number(roi_coding[roi_coding.count()-1].end()) + "' AND (gnomad IS NULL OR gnomad<=" + af + ") ORDER BY start");
			while(query.next())
			{
				QList<VariantTranscript> trans_infos;
				try
				{
					trans_infos = Variant::parseTranscriptString(query.value(5).toByteArray(), true);
				}
				catch(...) {} //do nothing (old RefSeq entries)

				foreach(const VariantTranscript& trans, trans_infos)
				{
					if ((trans.impact=="HIGH" || trans.impact=="MODERATE") &&  trans.id.startsWith(lct.name())) //no direct comparision of transcript name because we have mix transcripts with/without version number in NGSD.
					{
						QString variant_id = query.value(0).toString();
						QString var = chr.strNormalized(true) + ":" + query.value(1).toString() + "-" + query.value(2).toString() + " " + query.value(3).toString() + ">" + query.value(4).toString();
						QString af = query.value(6).toString();
						SqlQuery query2 = db.getQuery();
						query2.exec("SELECT s.disease_group, s.disease_status, s.id FROM sample s, processed_sample ps, project p, detected_variant dv, processing_system sys WHERE ps.processing_system_id=sys.id AND dv.processed_sample_id=ps.id AND ps.sample_id=s.id AND ps.project_id=p.id AND dv.variant_id=" + variant_id + " AND p.type='diagnostic' AND ps.quality!='bad' AND (sys.type='WES' OR sys.type='WGS')");

						QSet<int> sample_ids_done;
						int c_affected = 0;
						int c_unaffected = 0;
						int c_unknown = 0;
						QMap<QString, int> dg_affected;
						int samples_with_hpo = 0;
						QMap<QString, int> hpo_affected;
						while(query2.next())
						{
							QString disease_group = query2.value(0).toString();
							QString disease_status = query2.value(1).toString();
							int sample_id = query2.value(2).toInt();

							//skip duplicate samples and related samples
							if (sample_ids_done.contains(sample_id)) continue;
							sample_ids_done << sample_id;
							sample_ids_done.unite(db.relatedSamples(sample_id));

							if (disease_status=="Affected")
							{
								++c_affected;

								if (!dg_affected.contains(disease_group)) dg_affected[disease_group] = 0;
								dg_affected[disease_group] += 1;

								auto phenos = db.samplePhenotypes(QString::number(sample_id));
								if (phenos.count()>0) ++samples_with_hpo;
								foreach(const Phenotype& pheno, phenos)
								{
									QString hpo_name = pheno.name();
									if (!hpo_affected.contains(hpo_name)) hpo_affected[hpo_name] = 0;
									hpo_affected[hpo_name] += 1;
								}
							}
							else if (disease_status=="Unaffected")
							{
								++c_unaffected;
							}
							else
							{
								++c_unknown;
							}
						}
						if (c_affected<3) continue; //RESTRICTION at least 3 times in affected diagnostic WES/WGS samples
						QString dg_shared_by_affected;
						if (dg_affected.count()==1)
						{
							dg_shared_by_affected = dg_affected.keys().first();
						}
						QStringList hpos_shared_by_affected;
						foreach(QString hpo_name, hpo_affected.keys())
						{
							if (samples_with_hpo>=3 && hpo_affected[hpo_name]==samples_with_hpo) hpos_shared_by_affected << (hpo_name + " (" + QString::number(hpo_affected[hpo_name])+")");
							if (samples_with_hpo>=4 && hpo_affected[hpo_name]==samples_with_hpo-1) hpos_shared_by_affected << (hpo_name + " (" + QString::number(hpo_affected[hpo_name])+")");
						}
						stream << gene << "\t" << lct.name() << "\t" << var << "\t" << trans.hgvs_p << "\t" << trans.type << "\t" << trans.impact << "\t" << af << "\t" << db.getValue("SELECT class FROM variant_classification WHERE variant_id="+variant_id).toString() << "\t" << QString::number(c_affected) << "\t" << QString::number(c_unaffected) << "\t" << QString::number(c_unknown) << "\t" << dg_shared_by_affected << "\t" << QString::number(samples_with_hpo) << "\t" << hpos_shared_by_affected.join(", ") << "\n";
						stream.flush();
					}
				}
			}
		}
		*/

		//evaluation GSvar score/rank
		//init parameters and output stream
		bool test_domiant = false;
		bool test_recessive_hom = false; //in case of recessive - switch beteen hom and comp-het
		bool test_v1 = false;
		bool test_with_ngsd = false;
		bool test_with_clinvar = true;
		QString algorithm;
		QString suffix;
		if (test_v1)
		{
			algorithm = "GSvar_v1";
			suffix = test_domiant ? "_dominant" : "_recessive";
		}
		else if (test_domiant)
		{
			algorithm = "GSvar_v2_dominant";
		}
		else
		{
			algorithm = "GSvar_v2_recessive";
			suffix += test_recessive_hom ? "_hom" : "_comphet";
		}
		VariantScores::Parameters parameters;
		if (!test_with_ngsd)
		{
			parameters.use_ngsd_classifications = false;
			suffix += "_noNGSD";
		}
		if (!test_with_clinvar)
		{
			parameters.use_clinvar = false;
			suffix += "_noClinVar";
		}
		qDebug() << "algorithm:" << algorithm << "suffix:" << suffix;

		QSharedPointer<QFile> out_file = Helper::openFileForWriting("C:\\Marc\\ranking_" + QDate::currentDate().toString("yyyy-MM-dd") + "_" + algorithm + suffix + ".tsv");
		QTextStream out_stream(out_file.data());
		QStringList headers;
		headers << "ps" << "system" << "HPO" << "variant" << "variant_type"  << "filter" << "impact" << "inheritance" << "oe_lof" << "gnomAD" << "gnomAD sub" << "NGSD hom" << "NGSD het" << "NGSD mosaic" << "variants_scored" << "hpo_genes" << "score" << "score_explanations" << "rank";
		out_stream << "#" << headers.join("\t") << endl;

		QSet<QString> samples_skipped_other_variants;
		QSet<QString> samples_skipped_incorrect_causal_variant_count;
		QSet<QString> samples_skipped_no_hpo;
		int c_all = 0;
		int c_top1 = 0;
		int c_top3 = 0;
		int c_top10 = 0;
		int c_none = 0;
		QHash<QByteArray, GeneSet> pheno2genes_cache_;
		NGSD db;
		QString inheritance_constraint = test_domiant ? "AND (rcv.inheritance='AD' OR rcv.inheritance='XLD')" : "AND (rcv.inheritance='AR' OR rcv.inheritance='XLR')";
		QStringList ps_ids = db.getValues("SELECT DISTINCT ps.id FROM sample s, processed_sample ps, diag_status ds, report_configuration rc, report_configuration_variant rcv, project p, processing_system sys WHERE ps.processing_system_id=sys.id AND (sys.type='WGS' OR sys.type='WES') AND ps.project_id=p.id AND p.type='diagnostic' AND ps.sample_id=s.id AND ps.quality!='bad' AND ds.processed_sample_id=ps.id AND ds.outcome='significant findings' AND rc.processed_sample_id=ps.id AND rcv.report_configuration_id=rc.id AND rcv.causal='1' AND rcv.mosaic='0' AND rcv.type='diagnostic variant' AND s.disease_status='Affected' " + inheritance_constraint + " ORDER BY ps.id ASC");
		int ps_nr = 0;
		foreach(const QString& ps_id, ps_ids)
		{
			try
			{

				QString ps = db.processedSampleName(ps_id);
				int rc_id = db.reportConfigId(ps_id);

				//check that there are only small variants that are causal
				QList<int> causal_cnvs = db.getValuesInt("SELECT id FROM report_configuration_cnv WHERE report_configuration_id=:0 AND causal='1'", QString::number(rc_id));
				QList<int> causal_svs = db.getValuesInt("SELECT id FROM report_configuration_sv WHERE report_configuration_id=:0 AND causal='1'", QString::number(rc_id));
				QList<int> causal_other = db.getValuesInt("SELECT id FROM report_configuration_other_causal_variant WHERE report_configuration_id=:0", QString::number(rc_id));
				if (causal_cnvs.count()>0 || causal_svs.count()>0 || causal_other.count()>0)
				{
					samples_skipped_other_variants << ps;
					continue;
				}

				//check that there are the right number of causal variant for the inheritance mode (class 4/5, correct genotype, not mito)
				QString expected_inheritance = "(rcv.inheritance='AD' OR rcv.inheritance='XLD')";
				QString expected_genotype = "het";
				int expected_causal_variant_count = 1;
				if (!test_domiant && test_recessive_hom)
				{
					expected_inheritance = "(rcv.inheritance='AR' OR rcv.inheritance='XLR')";
					expected_genotype = "hom";
					expected_causal_variant_count = 1;
				}
				if (!test_domiant && !test_recessive_hom)
				{
					expected_inheritance = "(rcv.inheritance='AR' OR rcv.inheritance='XLR')";
					expected_genotype = "het";
					expected_causal_variant_count = 2;
				}
				VariantList causal_variants;
				foreach(QString v_id, db.getValues("SELECT rcv.variant_id FROM report_configuration_variant rcv, variant v, variant_classification vc WHERE rcv.variant_id=v.id AND rcv.variant_id=vc.variant_id AND rcv.report_configuration_id='" + QString::number(rc_id) + "' AND rcv.causal='1' AND rcv.type='diagnostic variant' AND v.chr!='chrMT' AND (vc.class='4' OR vc.class='5') AND " + expected_inheritance))
				{
					if (db.getValue("SELECT genotype FROM detected_variant WHERE processed_sample_id=" + ps_id + " AND variant_id=" + v_id + " AND mosaic=0").toByteArray()==expected_genotype)
					{
						causal_variants.append(db.variant(v_id));
					}
				}
				if (causal_variants.count()!=expected_causal_variant_count)
				{
					samples_skipped_incorrect_causal_variant_count << ps;
					continue;
				}

				//create phenotype list
				QHash<Phenotype, BedFile> phenotype_rois;
				QString sample_id = db.sampleId(ps);
				PhenotypeList phenotypes = db.getSampleData(sample_id).phenotypes;
				if (phenotypes.count()==0)
				{
					samples_skipped_no_hpo << ps;
					continue;
				}

				qDebug() << (++ps_nr) << ps;

				//determine HPO terms and HPO genes for output
				GeneSet hpo_genes;
				QStringList hpo_terms;
				foreach(const Phenotype& pheno, phenotypes)
				{
					//pheno > genes
					QByteArray pheno_accession = pheno.accession();
					if (!pheno2genes_cache_.contains(pheno_accession))
					{
						pheno2genes_cache_[pheno_accession] = db.phenotypeToGenes(db.phenotypeIdByAccession(pheno_accession), true);
					}
					GeneSet genes = pheno2genes_cache_[pheno_accession];

					hpo_terms << pheno.name();
					hpo_genes << genes;

					//genes > roi
					BedFile roi;
					foreach(const QByteArray& gene, genes)
					{
						if (!gene2region_cache_.contains(gene))
						{
							BedFile tmp = db.geneToRegions(gene, Transcript::ENSEMBL, "gene", true);
							tmp.clearAnnotations();
							tmp.merge();
							gene2region_cache_[gene] = tmp;
						}
						roi.add(gene2region_cache_[gene]);
					}
					roi.merge();

					phenotype_rois[pheno] = roi;
				}
				hpo_terms.sort();
				hpo_terms.removeDuplicates();

				//load variants
				QString gsvar = db.processedSamplePath(ps_id, PathType::GSVAR);
				QFileInfo gsvar_info(gsvar);
				if (gsvar_info.lastModified()<QDateTime::fromString("03.04.2023 14:10:00", "dd.MM.yyyy hh:mm:ss"))
				{
					qDebug() << "  skipped - GSvar file too old";
					continue;
				}

				//store local copy that is pre-filtered to reduce loading times
				QString tmp = "C:\\Marc\\ranking_prefiltered\\" + gsvar_info.fileName();
				if (!QFile::exists(tmp))
				{
					qDebug() << "  creating pre-filtered GSvar file...";
					VariantList variants_tmp;
					variants_tmp.load(gsvar);

					FilterCascade cascade = FilterCascade::fromText(QStringList() <<
						"Allele frequency	max_af=1.0" <<
						"Count NGSD	max_count=50	ignore_genotype=false	mosaic_as_het=false" <<
						"Annotated pathogenic	action=KEEP	sources=HGMD,ClinVar	also_likely_pathogenic=true" <<
						"Classification NGSD	action=KEEP	classes=4,5");
					FilterResult cascade_result = cascade.apply(variants_tmp);
					cascade_result.removeFlagged(variants_tmp);

					variants_tmp.store(tmp);
				}
				VariantList variants;
				variants.load(tmp);

				//score variants
				VariantScores::Result result = VariantScores::score(algorithm, variants, phenotype_rois, parameters);

				//prepare lambda for output
				int c_scored = VariantScores::annotate(variants, result, true);
				int i_filter = variants.annotationIndexByName("filter");
				int i_coding = variants.annotationIndexByName("coding_and_splicing");
				int i_gnomad = variants.annotationIndexByName("gnomAD");
				int i_gnomad_sub = variants.annotationIndexByName("gnomAD_sub");
				int i_rank = variants.annotationIndexByName("GSvar_rank");
				int i_score = variants.annotationIndexByName("GSvar_score");
				int i_score_exp = variants.annotationIndexByName("GSvar_score_explanations");

				auto addLine = [&] (const Variant& v, QString variant_type)
				{
					QList<int> affected_cols = variants.getSampleHeader().sampleColumns(true);
					if (affected_cols.count()!=1) THROW(ArgumentException, "VariantScores: Algorihtm 'GSvar_v1' can only be applied to variant lists with exactly one affected patient!");
					int i_genotye = affected_cols[0];
					QString genotype = v.annotations()[i_genotye].trimmed();

					QString gnomad_af = v.annotations()[i_gnomad].trimmed();
					if (gnomad_af.isEmpty()) gnomad_af = "0";

					QString gnomad_af_sub = v.annotations()[i_gnomad_sub].trimmed();
					if (gnomad_af_sub.isEmpty())
					{
						gnomad_af_sub = "0";
					}
					else
					{
						double max = 0;
						QStringList parts = gnomad_af_sub.split(",");
						foreach(const QString& part, parts)
						{
							double value = Helper::toDouble(part, "gnomAD sub-population AF");
							if (value>max) max = value;
						}
						gnomad_af_sub = QString::number(max,'f', 5);
					}

					QString rank = v.annotations()[i_rank].trimmed();
					if (rank.isEmpty()) rank = "-3";

					QString score = v.annotations()[i_score].trimmed();
					if (score.isEmpty()) rank = "-3";

					QByteArrayList impact;
					QByteArrayList gene_inheritance;
					QByteArrayList oe_lof;
					QMap<QByteArray, QSet<QByteArray>> gene2impact;
					QList<VariantTranscript> transcript_info = v.transcriptAnnotations(i_coding);
					foreach(const VariantTranscript& transcript, transcript_info)
					{
						gene2impact[transcript.gene] << transcript.impact;
					}
					for (auto it=gene2impact.begin(); it!=gene2impact.end(); ++it)
					{
						QByteArray gene = it.key();

						impact += gene + "=" + it.value().toList().join("/");

						GeneInfo gene_info = db.geneInfo(gene);
						gene_inheritance += gene + "=" + gene_info.inheritance.toLatin1();
						oe_lof += gene + "=" + gene_info.oe_lof.toLatin1();
					}

					//get system name
					int sys_id = db.processingSystemIdFromProcessedSample(ps);
					QString system = db.getProcessingSystemData(sys_id).name;

					//get variant counts
					QString var_id = db.variantId(v);
					GenotypeCounts geno_counts = db.genotypeCountsCached(var_id);

					QStringList cols;
					cols << ps << system << hpo_terms.join(", ") << (v.toString() + " (" + genotype + ")") << variant_type  << v.annotations()[i_filter]<< impact.join(", ") << gene_inheritance.join(", ") << oe_lof.join(", ") << gnomad_af << gnomad_af_sub << QString::number(geno_counts.hom) << QString::number(geno_counts.het) << QString::number(geno_counts.mosaic) << QString::number(c_scored) << QString::number(hpo_genes.count()) << score << v.annotations()[i_score_exp] << rank;
					out_stream << cols.join("\t") << endl;
				};

				//determine causal variant ranks to fix statistics for comp-het case
				QList<int> causal_variant_ranks;
				if (!test_domiant && !test_recessive_hom)
				{
					for (int i=0; i<causal_variants.count(); ++i)
					{

						try
						{

						int var_index = variants.indexOf(causal_variants[i]);
						if (var_index==-1) THROW(Exception, "not found");

						int rank = Helper::toInt(variants[var_index].annotations()[i_rank]);
						causal_variant_ranks << rank;
						}
						catch(Exception& e)
						{
							causal_variant_ranks << 9999;
						}
					}
				}
				std::sort(causal_variant_ranks.begin(), causal_variant_ranks.end());

				//perform rank statistics on causal variants
				QSet<int> noncausal_indices_added;
				for (int i=0; i<causal_variants.count(); ++i)
				{
					++c_all;
					Variant causal_var = causal_variants[i];
					//qDebug() << __LINE__ << ps << i << causal_var.toString();

					int var_index = variants.indexOf(causal_var);
					if (var_index==-1)
					{
						causal_var.annotations() << QVector<QByteArray>(100).toList(); //no annotations would crash
						addLine(causal_var, "missed causal variant - not in pre-filtered variant list");
						continue;
					}

					const Variant& var = variants[var_index];

					addLine(var, "causal " + QString::number(i+1) + "/" + QString::number(causal_variants.count()));

					//statistics
					try
					{
						int rank = Helper::toInt(var.annotations()[i_rank]);
						if (!test_domiant && !test_recessive_hom && rank>causal_variant_ranks[0]) rank -= 1; //ignore first variant in comp-het case to make the comparison between inheritance modes fair
						if (rank==1) ++c_top1;
						if (rank<=3) ++c_top3;
						if (rank<=10) ++c_top10;

						//store top 5 variants ranking higher than the causal variant
						if (rank>1 && rank<=10)
						{
							for(int v=0; v<variants.count() ; ++v)
							{
								const Variant& var2 = variants[v];
								if (causal_variants.contains(var2)) continue;
								if (noncausal_indices_added.contains(v)) continue;

								bool ok = false;
								int rank2 = variants[v].annotations()[i_rank].toInt(&ok);
								if (ok && rank2<rank && rank2<=5)
								{
									addLine(variants[v], "non-causal high-ranked variant");
									noncausal_indices_added << v;
								}
							}
						}
					}
					catch(Exception& e)
					{
						addLine(var, "missed causal variant - not ranked");
						++c_none;
					}
				}
			}
			catch(Exception& e)
			{
				qDebug() << "  Error processing sample:" << e.message();
				continue;
			}
		}
		out_stream << "##Number of samples skipped because of other causal variant: " << samples_skipped_other_variants.count() << endl;
		out_stream << "##Number of samples skipped because of incorrect causal variant count: " << samples_skipped_incorrect_causal_variant_count.count() << endl;
		out_stream << "##Number of samples skipped because of missing HPO terms: " << samples_skipped_no_hpo.count() << endl;
		out_stream << "##Number of samples used for benchmark: " << ps_nr << endl;
		out_stream << "##Number of variants used for benchmark: " << c_all << endl;
		out_stream << "##Rank1: " << QString::number(c_top1) << " (" + QString::number(100.0*c_top1/c_all, 'f', 2) << "%)" << endl;
		out_stream << "##Top3 : " << QString::number(c_top3) << " (" + QString::number(100.0*c_top3/c_all, 'f', 2) << "%)" << endl;
		out_stream << "##Top10: " << QString::number(c_top10) << " (" + QString::number(100.0*c_top10/c_all, 'f', 2) << "%)" << endl;
		out_stream << "##None : " << QString::number(c_none) << " (" + QString::number(100.0*c_none/c_all, 'f', 2) << "%)" << endl;

		//import of sample relations from GenLab
		/*
		QStringList pairs;
		pairs << "DX070696	DX070760";

		NGSD db;
		GenLabDB db_genlab;
		foreach(QString pair, pairs)
		{
			QStringList parts = pair.split("\t");
			if (parts.count()!=2)
			{
				qDebug() << "Error: invalid line: " << pair;
				break;
			}
			QString sample1 = parts[0];
			QString sample2 = parts[1];
			qDebug() << sample1 << sample2;

			//check one direction
			QList<SampleRelation> relatives = db_genlab.relatives(sample1);
			foreach(const SampleRelation& rel, relatives)
			{
				if (rel.sample1==sample2)
				{
					db.addSampleRelation(rel);
				}
			}

			//check other direction
			relatives = db_genlab.relatives(sample2);
			foreach(SampleRelation rel, relatives)
			{
				if (rel.sample1==sample1)
				{
					db.addSampleRelation(rel);
				}
			}
		}
		*/

		//non-causal variants annotation
		/*
		NGSD db;
		QStringList input;
		input << TODO
		foreach(QString ps, input)
		{
			QString ps_id = db.processedSampleId(ps);
			QString text;
			int rc_id = db.reportConfigId(ps_id);
			if (rc_id!=-1)
			{
				//find causal small variants
				QStringList causal_ids = db.getValues("SELECT variant_id FROM report_configuration_variant WHERE causal='0' AND exclude_artefact='0' AND exclude_frequency='0' AND exclude_phenotype='0' AND exclude_mechanism='0' AND exclude_other='0' AND report_configuration_id=" + QString::number(rc_id));
				foreach(QString id, causal_ids)
				{
					Variant var = db.variant(id);
					QString genotype = db.getValue("SELECT genotype FROM detected_variant WHERE processed_sample_id='" + ps_id + "' AND variant_id='" + id + "'").toString();
					QString genes = db.genesOverlapping(var.chr(), var.start(), var.end(), 5000).join(", ");
					QString var_class = db.getValue("SELECT class FROM variant_classification WHERE variant_id='" + id + "'").toString();
					text += ", small variant: " + var.toString() + " (genotype:" + genotype + " genes:" + genes;
					if (var_class != "") text += " classification:" + var_class; // add classification, if exists
					text += ")";
				}

				//find causal CNVs
				causal_ids = db.getValues("SELECT cnv_id FROM report_configuration_cnv WHERE  causal='0' AND exclude_artefact='0' AND exclude_frequency='0' AND exclude_phenotype='0' AND exclude_mechanism='0' AND exclude_other='0' AND report_configuration_id=" + QString::number(rc_id));
				foreach(QString id, causal_ids)
				{
					CopyNumberVariant var = db.cnv(id.toInt());
					QString cn = db.getValue("SELECT cn FROM cnv WHERE id='" + id + "'").toString();
					QString cnv_class = db.getValue("SELECT class FROM report_configuration_cnv WHERE cnv_id='" + id + "'", false).toString();
					text += ", CNV: " + var.toString() + " (cn:" + cn;
					if (cnv_class != "") text += " classification:" + cnv_class; // add classification, if exists
					text += ")";
				}

				//find causal SVs
				QStringList sv_id_columns = QStringList() << "sv_deletion_id" << "sv_duplication_id" << "sv_insertion_id" << "sv_inversion_id" << "sv_translocation_id";
				QList<StructuralVariantType> sv_types = {StructuralVariantType::DEL, StructuralVariantType::DUP, StructuralVariantType::INS, StructuralVariantType::INV, StructuralVariantType::BND};
				BedpeFile svs;
				for (int i = 0; i < sv_id_columns.size(); ++i)
				{
					causal_ids = db.getValues("SELECT " + sv_id_columns.at(i) + " FROM report_configuration_sv WHERE  causal='0' AND exclude_artefact='0' AND exclude_frequency='0' AND exclude_phenotype='0' AND exclude_mechanism='0' AND exclude_other='0' AND report_configuration_id=" + QString::number(rc_id) + " AND " + sv_id_columns.at(i) + " IS NOT NULL");

					foreach(QString id, causal_ids)
					{
						BedpeLine var = db.structuralVariant(id.toInt(), sv_types.at(i), svs, true);
						QString sv_class = db.getValue("SELECT class FROM report_configuration_sv WHERE " + sv_id_columns[i] + "='" + id + "'", false).toString();
						text += ", SV: " + var.toString();
						if (sv_class != "") text += " (classification:" + sv_class + ")"; // add classification, if exists
					}
				}
			}
			qDebug() << ps << "\t" << text;
		}
		*/

		//Export GenLab dates for reanalysis of unsolved samples
		/*
		TsvFile output;
		output.addHeader("ps");
		output.addHeader("yearOfBirth");
		output.addHeader("yearOfOrderEntry");
		GenLabDB db;
		TsvFile file;
		file.load("W:\\share\\evaluations\\2020_07_14_reanalysis_pediatric_cases\\samples.tsv");
		int i=0;
		QStringList ps_names = file.extractColumn(1);
		foreach(QString ps, ps_names)
		{
			qDebug() << ++i << "/" << ps_names.count() << ps;
			output.addRow(QStringList() << ps << db.yearOfBirth(ps) << db.yearOfOrderEntry(ps));
		}
		output.store("W:\\share\\evaluations\\2020_07_14_reanalysis_pediatric_cases\\+documentation\\genlab_export_dates_" + QDate::currentDate().toString("yyyy_MM_dd")+".tsv");
		*/

		//import preferred transcripts
		/*
		NGSD db;
		QString filename = GSvarHelper::applicationBaseName() + "_preferred_transcripts.tsv";
		QStringList lines = Helper::loadTextFile(filename, true, '#', true);
		foreach(const QString& line, lines)
		{
			QByteArrayList parts = line.toUtf8().replace(',', '\t').split('\t');
			if (parts.count()>=2)
			{
				QByteArray gene = parts[0].trimmed();
				for (int i=1; i<parts.count(); ++i)
				{
					QByteArray transcript = parts[i].trimmed();
					qDebug() << gene << transcript;
					try
					{
						qDebug() << "  success: " << db.addPreferredTranscript(transcript);
					}
					catch(Exception& e)
					{

						qDebug() << "  failed: " << e.message();
					}
				}
			}
		}
		*/

		//import sample meta data from GenLab
		/*
		GenLabDB genlab;
		NGSD db;
		ProcessedSampleSearchParameters params;
		params.p_type = "diagnostic";
		params.sys_type = "WGS";
		params.include_bad_quality_samples = false;
		params.include_tumor_samples = false;
		params.include_merged_samples = false;
		params.include_bad_quality_runs = false;
		params.run_finished = true;
		DBTable ps_table = db.processedSampleSearch(params);
		QStringList ps_list = ps_table.extractColumn(0);

		int ps_start_index = -1;
		int i=0;
		foreach(QString ps, ps_list)
		{
			++i;
			if (i<ps_start_index) continue;

			qDebug() << i << "/" << ps_list.size() << " - " << ps;
			genlab.addMissingMetaDataToNGSD(ps, true, true, true, true, false);
		}
		*/

		//initial import of patient identifiers from GenLab (diagnostic samples only)
		/*
		NGSD db;
		GenLabDB db_genlab;
		SqlQuery query = db.getQuery();
		query.exec("SELECT s.id, concat(s.name, '_0', ps.process_id), s.patient_identifier FROM sample s, processed_sample ps, project p WHERE s.id=ps.sample_id AND p.id=ps.project_id AND p.type='diagnostic' ORDER BY ps.id ASC");
		while(query.next())
		{
			QString s_id = query.value(0).toString().trimmed();
			QString ps = query.value(1).toString().trimmed();
			QString patient_id_old = query.value(2).toString().trimmed();

			QString patient_id = db_genlab.patientIdentifier(ps);
			if (patient_id=="") continue;

			//check for mismatches
			if (patient_id_old!="")
			{
				if (patient_id!=patient_id_old) qDebug() << "MISMATCH:" << ps << "NGSD=" << patient_id_old << "GenLab=" << patient_id;
				continue;
			}

			qDebug() << "UPDATE:" << ps << patient_id;
			db.getQuery().exec("UPDATE sample SET patient_identifier='" + patient_id + "' WHERE id='" + s_id + "'");
		}
		*/

		//replace obsolete terms used in disease_info table
		/*
		OntologyTermCollection terms("W:\\GRCh38\\share\\data\\dbs\\HPO\\hp.obo", false);
		int c_obsolote = 0;
		int c_replaced = 0;
		QByteArrayList invalid;
		NGSD db;
		SqlQuery query = db.getQuery();
		query.exec("SELECT id, disease_info FROM sample_disease_info WHERE type='HPO term id' order by disease_info ASC");
		while(query.next())
		{
			QByteArray hpo_id = query.value("disease_info").toByteArray().trimmed();
			try
			{
				const OntologyTerm& term = terms.getByID(hpo_id);
				if (term.isObsolete())
				{
					++c_obsolote;
					if (!term.replacedById().isEmpty())
					{
						qDebug() << term.id() << " > " << term.replacedById();
						db.getQuery().exec("UPDATE sample_disease_info SET disease_info='" + term.replacedById() + "' WHERE id=" + query.value("id").toByteArray());
						++c_replaced;
					}
				}
			}
			catch(const ArgumentException&)
			{
				if (!invalid.contains(hpo_id)) invalid << hpo_id;
			}
		}
		if (c_obsolote>0)
		{
			qDebug() << "Found " << c_obsolote << " obsolete HPO terms in table 'disease_info'. Replaced " << c_replaced << " of these!" << endl;
		}
		if (invalid.count()>0)
		{
			qDebug() << "Found " << invalid.count() << " invalid HPO terms in table 'disease_info': '" << invalid.join("', '") << "'" << endl;
		}
		*/

		//search for SVs with breakpoints inside genes matching the phenotype
		/*
		NGSD db;
		QSharedPointer<QFile> file = Helper::openFileForWriting("C:\\Marc\\large_sv_breakpoints.tsv");
		QTextStream ostream(file.data());

		QSharedPointer<QFile> file2 = Helper::openFileForWriting("C:\\Marc\\large_sv_breakpoints_stats.tsv");
		QTextStream ostream2(file2.data());

		QStringList ps_ids = db.getValues("SELECT ps.id FROM processed_sample ps, processing_system sys, project p, sample s WHERE ps.processing_system_id=sys.id AND sys.name_short='TruSeqPCRfree' AND ps.project_id=p.id AND p.type='diagnostic' AND ps.sample_id=s.id AND s.tumor='0' AND ps.quality!='bad' AND ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples)");
		ostream2 << "##diagnostic germline genomes: " << ps_ids.count() << endl;

		QMap<QString, QStringList> ps_id2hpos;
		foreach(QString ps_id, ps_ids)
		{
			QStringList hpo_terms = db.getValues("SELECT disease_info FROM sample_disease_info sdi, processed_sample ps WHERE sdi.sample_id=ps.sample_id AND ps.id='"+ps_id+"' AND type='HPO term id'");
			std::for_each(hpo_terms.begin(), hpo_terms.end(), [](QString& value){ value = value.toUpper().trimmed(); });
			hpo_terms.removeAll("HP:0032322"); //healthy
			if (hpo_terms.count()>0)
			{
				ps_id2hpos[ps_id] = hpo_terms;
			}
		}
		ostream2 << "##diagnostic germline genomes with HPO terms: " << ps_id2hpos.count() << endl;

		FilterCascade filters = FilterCascade::fromText(QStringList()   << "SV type	Structural variant type=INV"
																		<< "SV remove chr type	chromosome type=special chromosomes"
																		<< "SV allele frequency NGSD	max_af=0.1"
																		<< "SV filter columns	entries=PASS	action=FILTER"
																		<< "SV break point density NGSD	max_density=20	remove_strict=no"
																		<< "SV PE read depth	PE Read Depth=5	only_affected=no"
																		<< "SV size	min_size=500000	max_size=0");
		ostream << "#ps" << "\t" << "variant" << "\t" << "report_config_of_variant" << "\t" << "disease" << "\t" << "outcome" << "\t" << "user_rc_create" << "\t" << "user_rc_last_edit" << endl;
		ostream2 << "#ps" << "\t" << "HPO_terms" << "\t" << "HPO_genes" << "\t" << "HPO_roi_bases" << "\t" << "SVs" << "\t" << "SVs_after_filter" << "\t" << "SVs_in_roi" << endl;
		for(auto it=ps_id2hpos.begin(); it!=ps_id2hpos.end(); ++it)
		{
			QString ps_id = it.key();
			QStringList hpos = it.value();
			QString ps = db.processedSampleName(ps_id);
			int rc_id = db.reportConfigId(ps_id);
			QString disease = db.getSampleData(db.sampleId(ps)).disease_group + " (" + db.getSampleData(db.sampleId(ps)).disease_status + ")";
			QString outcome = db.getDiagnosticStatus(ps_id).outcome;
			QString user_create = db.getValue("SELECT u.name FROM report_configuration rc, user u WHERE rc.created_by=u.id AND rc.id=" + QString::number(rc_id)).toString();
			QString user_last = db.getValue("SELECT u.name FROM report_configuration rc, user u WHERE rc.last_edit_by=u.id AND rc.id=" + QString::number(rc_id)).toString();

			GeneSet pheno_genes;
			foreach(const QString& hpo, hpos)
			{
				int pheno_id = db.phenotypeIdByAccession(hpo.toUtf8(), false);
				if (pheno_id==-1) continue;
				Phenotype pheno = db.phenotype(pheno_id);
				GeneSet genes = db.phenotypeToGenesbySourceAndEvidence(db.phenotypeIdByAccession(pheno.accession()), QSet<PhenotypeSource>(), QSet<PhenotypeEvidenceLevel>(), true);
				pheno_genes << genes;
			}

			//convert genes to ROI (using a cache to speed up repeating queries)
			BedFile roi;
			foreach(const QByteArray& gene, pheno_genes)
			{
				if (!gene2region_cache_.contains(gene))
				{
					BedFile tmp = db.geneToRegions(gene, Transcript::ENSEMBL, "gene", true);
					tmp.clearAnnotations();
					tmp.merge();
					gene2region_cache_[gene] = tmp;
				}
				roi.add(gene2region_cache_[gene]);
			}
			roi.merge();
			if (roi.baseCount()==0) continue;

			QVariant callset_id = db.getValue("SELECT id FROM sv_callset WHERE processed_sample_id='" + ps_id + "'", true);
			if (!callset_id.isValid()) continue;

			//get rare SVs of this case
			QString sv_file = db.processedSamplePath(ps_id, PathType::STRUCTURAL_VARIANTS);
			if (!QFile::exists(sv_file)) continue;

			BedpeFile svs;
			svs.load(sv_file);
			int c_svs_nofilter = svs.count();

			FilterResult filter_res = filters.apply(svs, false, false);
			filter_res.removeFlagged(svs);

			int c_snv_in_roi = 0;
			for(int i=0; i<svs.count(); ++i)
			{
				const BedpeLine& sv = svs[i];

				bool overlap = roi.overlapsWith(sv.chr1(), sv.start1(), sv.end1()) || roi.overlapsWith(sv.chr2(), sv.start2(), sv.end2());
				if (overlap)
				{
					//check if RC exists
					QString rc_exists = "no";

					if (rc_id!=-1)
					{
						QString sv_id = db.svId(sv, callset_id.toInt(), svs);

						QVariant causal = db.getValue("SELECT causal FROM report_configuration_sv WHERE report_configuration_id='" + QString::number(rc_id) + "' AND sv_inversion_id='" + sv_id + "'", true);
						if (causal.isValid())
						{
							rc_exists = "yes";
							if (causal.toBool())
							{
								rc_exists += " - causal";
							}
						}
					}

					ostream << ps << "\t" << sv.toString() << "\t" << rc_exists << "\t" << disease << "\t" << outcome << "\t" << user_create << "\t" << user_last << endl;
					c_snv_in_roi += 1;

				}
			}
			ostream2 << ps << "\t" << hpos.join(", ") << "\t" << pheno_genes.count() << "\t" << roi.baseCount() << "\t" << c_svs_nofilter << "\t" << svs.count() << "\t" << c_snv_in_roi << endl;
		}
		*/
	}
	else if (user=="ahgscha1")
	{
	}
	else if (user=="ahschul1")
	{
		//Test sample sheet
		QString run_id = NGSD().getValue("SELECT id FROM sequencing_run WHERE name='#03178'").toString();

		QMessageBox::information(this, "SampleSheet", NGSD().createSampleSheet(run_id.toInt()));

	}
	else if (user=="ahott1a1")
	{

	}

	qDebug() << "Elapsed time debugging:" << Helper::elapsedTime(timer, true);
}



void MainWindow::on_actionConvertVcfToGSvar_triggered()
{
	VariantConversionWidget* widget = new VariantConversionWidget();
	widget->setMode(VariantConversionWidget::VCF_TO_GSVAR);
	auto dlg = GUIHelper::createDialog(widget, "Variant conversion (VCF > GSvar)");
	addModelessDialog(dlg);
}

void MainWindow::on_actionConvertHgvsToGSvar_triggered()
{
	VariantConversionWidget* widget = new VariantConversionWidget();
	widget->setMode(VariantConversionWidget::HGVSC_TO_GSVAR);
	auto dlg = GUIHelper::createDialog(widget, "Variant conversion (HGVS.c > GSvar)");
	addModelessDialog(dlg);
}

void MainWindow::on_actionConvertGSvarToVcf_triggered()
{
	VariantConversionWidget* widget = new VariantConversionWidget();
	widget->setMode(VariantConversionWidget::GSVAR_TO_VCF);
	auto dlg = GUIHelper::createDialog(widget, "Variant conversion (GSvar > VCF)");
	addModelessDialog(dlg);
}

void MainWindow::on_actionCytobandsToRegions_triggered()
{
	CytobandToRegionsDialog dlg(this);

	dlg.exec();
}

void MainWindow::on_actionRegionToGenes_triggered()
{
	QString title = "Region > Genes";

	try
	{
		//get region string
		QString region_text = QInputDialog::getText(this, title, "genomic region");
		if (region_text=="") return;

		QApplication::setOverrideCursor(Qt::BusyCursor);

		//convert to region
		Chromosome chr;
		int start, end;
		NGSHelper::parseRegion(region_text, chr, start, end);

		//convert region to gene set
		NGSD db;
		GeneSet genes = db.genesOverlapping(chr, start, end);

		QApplication::restoreOverrideCursor();

		//show results
		ScrollableTextDialog dlg(this, title);
		dlg.setReadOnly(true);
		dlg.setWordWrapMode(QTextOption::NoWrap);
		dlg.appendLine("#GENE\tOMIM_GENE\tOMIM_PHENOTYPES");
		foreach (const QByteArray& gene, genes)
		{
			QList<OmimInfo> omim_genes = db.omimInfo(gene);
			foreach (const OmimInfo& omim_gene, omim_genes)
			{
				dlg.appendLine(gene + "\t" + omim_gene.gene_symbol + "\t" + omim_gene.phenotypes.toString());
			}
		}
		dlg.exec();
	}
	catch(Exception& e)
	{
		QApplication::restoreOverrideCursor();
		QMessageBox::warning(this, title, "Error:\n" + e.message());
		return;
	}
}

void MainWindow::on_actionSearchSNVs_triggered()
{
	SmallVariantSearchWidget* widget = new SmallVariantSearchWidget();
	auto dlg = GUIHelper::createDialog(widget, "Small variants search");
	addModelessDialog(dlg);
}

void MainWindow::on_actionSearchCNVs_triggered()
{
	CnvSearchWidget* widget = new CnvSearchWidget();
	auto dlg = GUIHelper::createDialog(widget, "CNV search");
	addModelessDialog(dlg);
}

void MainWindow::on_actionSearchSVs_triggered()
{
	SvSearchWidget* widget = new SvSearchWidget();
	auto dlg = GUIHelper::createDialog(widget, "SV search");
	addModelessDialog(dlg);
}

void MainWindow::on_actionUploadVariantToClinVar_triggered()
{
	QSharedPointer<QDialog> dlg = QSharedPointer<QDialog>(new ClinvarUploadDialog(this));
	addModelessDialog(dlg);
}

void MainWindow::on_actionShowPublishedVariants_triggered()
{
	PublishedVariantsWidget* widget = new PublishedVariantsWidget();

	auto dlg = GUIHelper::createDialog(widget, "Published variants");
	dlg->exec();
}

void MainWindow::on_actionAlleleBalance_triggered()
{
	AlleleBalanceCalculator* widget = new AlleleBalanceCalculator();
	auto dlg = GUIHelper::createDialog(widget, "Allele balance of heterzygous variants");
	dlg->exec();
}

void MainWindow::on_actionLiftOver_triggered()
{
	LiftOverWidget* widget = new LiftOverWidget(this);
	auto dlg = GUIHelper::createDialog(widget, "Lift-over genome coordinates");
	addModelessDialog(dlg);
}

void MainWindow::on_actionGetGenomicSequence_triggered()
{
	QString title = "Get genomic sequence";
	try
	{
		//get region
		QString region_text = QInputDialog::getText(this, title, "genomic region:");
		if (region_text=="") return;

		Chromosome chr;
		int start, end;
		NGSHelper::parseRegion(region_text, chr, start, end);

		//get sequence
		QString genome_file = Settings::string("reference_genome", false);
		FastaFileIndex genome_idx(genome_file);
		int length = end-start+1;
		Sequence sequence = genome_idx.seq(chr, start, length, true);

		//copy to clipboard
		QApplication::clipboard()->setText(sequence);

		//show message
		if (sequence.length()>100)
		{
			sequence.resize(100);
			sequence += "...";
		}
		QMessageBox::information(this, title, "Extracted reference sequence of region " + chr.strNormalized(true) + ":" + QString::number(start) + "-" + QString::number(end) + " (length " + QString::number(length) + "):\n" + sequence + "\n\nThe sequence was copied to the clipboard.");
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, title, "Error getting reference sequence:\n" + e.message());
	}
}

void MainWindow::on_actionBlatSearch_triggered()
{
	BlatWidget* widget = new BlatWidget(this);

	auto dlg = GUIHelper::createDialog(widget, "BLAT search");
	addModelessDialog(dlg);
}


void MainWindow::on_actionClose_triggered()
{
	loadFile();
}

void MainWindow::on_actionCloseMetaDataTabs_triggered()
{
	for (int t=ui_.tabs->count()-1; t>0; --t)
	{
		closeTab(t);
	}
}

void MainWindow::on_actionIgvClear_triggered()
{
	IgvSessionManager::get(0).clear();
}

void MainWindow::on_actionIgvDocumentation_triggered()
{
	QDesktopServices::openUrl(QUrl("https://software.broadinstitute.org/software/igv/UserGuide"));
}

void MainWindow::on_actionSV_triggered()
{
	if(filename_ == "") return;

	if (!svs_.isValid())
	{
		QMessageBox::information(this, "SV file missing", "No structural variant file is present in the analysis folder!");
		return;
	}

	//create list of genes with heterozygous variant hits
	GeneSet het_hit_genes;
	int i_genes = variants_.annotationIndexByName("gene", true, false);
	QList<int> i_genotypes = variants_.getSampleHeader().sampleColumns(true);
	i_genotypes.removeAll(-1);

	if (i_genes!=-1 && i_genotypes.count()>0)
	{
		//check that a filter was applied (otherwise this can take forever)
		int passing_vars = filter_result_.countPassing();
		if (passing_vars>3000)
		{
			int res = QMessageBox::question(this, "Continue?", "There are " + QString::number(passing_vars) + " variants that pass the filters.\nGenerating the list of candidate genes for compound-heterozygous hits may take very long for this amount of variants.\nDo you want to continue?", QMessageBox::Yes, QMessageBox::No);
			if(res==QMessageBox::No) return;
		}
		for (int i=0; i<variants_.count(); ++i)
		{
			if (!filter_result_.passing(i)) continue;

			bool all_genos_het = true;
			foreach(int i_genotype, i_genotypes)
			{
				if (variants_[i].annotations()[i_genotype]!="het")
				{
					all_genos_het = false;
				}
			}
			if (!all_genos_het) continue;
			het_hit_genes.insert(GeneSet::createFromText(variants_[i].annotations()[i_genes], ','));
		}
	}
	else if (variants_.type()!=SOMATIC_PAIR && variants_.type() != SOMATIC_SINGLESAMPLE)
	{
		QMessageBox::information(this, "Invalid variant list", "Column for genes or genotypes not found in variant list. Cannot apply compound-heterozygous filter based on variants!");
	}

	try
	{
		//determine processed sample ID (needed for report config)
		QString ps_id = "";
		QSharedPointer<ReportConfiguration> report_config = nullptr;
		if (germlineReportSupported())
		{
			ps_id = NGSD().processedSampleId(germlineReportSample(), false);
			report_config = report_settings_.report_config;
		}

		//open SV widget
		SvWidget* list;
		if(svs_.isSomatic())
		{
			// somatic
			list = new SvWidget(svs_, ps_id, ui_.filters, het_hit_genes, gene2region_cache_, this);
		}
		else
		{
			// germline single, trio or multi sample
			list = new SvWidget(svs_, ps_id, ui_.filters, report_config, het_hit_genes, gene2region_cache_, this);
		}

		auto dlg = GUIHelper::createDialog(list, "Structural variants of " + variants_.analysisName());
		addModelessDialog(dlg);
	}
	catch(FileParseException error)
	{
		QMessageBox::warning(this,"File Parse Exception",error.message());
	}
	catch(FileAccessException error)
	{
		QMessageBox::warning(this,"SV file not found",error.message());
	}
}

void MainWindow::on_actionCNV_triggered()
{
	if (filename_=="") return;

	if (!cnvs_.isValid())
	{
		QMessageBox::information(this, "CNV file missing", "No copy-number file is present in the analysis folder!");
		return;
	}

	AnalysisType type = variants_.type();

	//create list of genes with heterozygous variant hits
	GeneSet het_hit_genes;
	int i_genes = variants_.annotationIndexByName("gene", true, false);
	QList<int> i_genotypes = variants_.getSampleHeader().sampleColumns(true);
	i_genotypes.removeAll(-1);

	if (i_genes!=-1 && i_genotypes.count()>0)
	{
		//check that a filter was applied (otherwise this can take forever)
		int passing_vars = filter_result_.countPassing();
		if (passing_vars>3000)
		{
			int res = QMessageBox::question(this, "Continue?", "There are " + QString::number(passing_vars) + " variants that pass the filters.\nGenerating the list of candidate genes for compound-heterozygous hits may take very long for this amount of variants.\nPlease set a filter for the variant list, e.g. the recessive filter, and retry!\nDo you want to continue?", QMessageBox::Yes, QMessageBox::No);
			if(res==QMessageBox::No) return;
		}
		for (int i=0; i<variants_.count(); ++i)
		{
			if (!filter_result_.passing(i)) continue;

			bool all_genos_het = true;
			foreach(int i_genotype, i_genotypes)
			{
				if (variants_[i].annotations()[i_genotype]!="het")
				{
					all_genos_het = false;
				}
			}
			if (!all_genos_het) continue;
			het_hit_genes.insert(GeneSet::createFromText(variants_[i].annotations()[i_genes], ','));
		}
	}
	else if (type!=SOMATIC_PAIR && type!=SOMATIC_SINGLESAMPLE)
	{
		QMessageBox::information(this, "Invalid variant list", "Column for genes or genotypes not found in variant list. Cannot apply compound-heterozygous filter based on variants!");
	}

	//determine processed sample ID (needed for report config)
	QString ps_id = "";
	if (germlineReportSupported())
	{
		ps_id = NGSD().processedSampleId(germlineReportSample(), false);
	}

	CnvWidget* list;
	if(cnvs_.type() == CnvListType::CLINCNV_TUMOR_NORMAL_PAIR || cnvs_.type() == CnvListType::CLINCNV_TUMOR_ONLY)
	{
		list = new CnvWidget(cnvs_, ps_id, ui_.filters, somatic_report_settings_.report_config, het_hit_genes, gene2region_cache_);
		connect(list, SIGNAL(storeSomaticReportConfiguration()), this, SLOT(storeSomaticReportConfig()));
	}
	else
	{
		list = new CnvWidget(cnvs_, ps_id, ui_.filters, report_settings_.report_config, het_hit_genes, gene2region_cache_);
	}

	auto dlg = GUIHelper::createDialog(list, "Copy number variants of " + variants_.analysisName());
	addModelessDialog(dlg);

	//mosaic CNVs
	if (type==GERMLINE_SINGLESAMPLE)
	{
		FileLocation mosaic_file = GlobalServiceProvider::fileLocationProvider().getAnalysisMosaicCnvFile();
		if (mosaic_file.exists)
		{
			QStringList mosaic_data = Helper::loadTextFile(mosaic_file.filename, false, '#', true);
			if (!mosaic_data.isEmpty())
			{
				ScrollableTextDialog dlg(this, "Possible mosaic CNV(s) detected!");
				dlg.appendLine("#CHR\tSTART\tEND\tCOPY NUMBER");

				foreach (const QString& line, mosaic_data)
				{
					if(line.trimmed().isEmpty() || line.startsWith("#")) continue;

					QStringList parts = line.split("\t");
					if(parts.length()<4)
					{
						Log::warn("Mosaic CNV file '" + mosaic_file.filename + "' has line with less than 4 elements: " + line);
					}
					else
					{
						dlg.appendLine(parts.mid(0, 4).join("\t"));
					}
				}
				dlg.exec();
			}
		}
	}
}

void MainWindow::on_actionROH_triggered()
{
	if (filename_=="") return;

	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE) return;

	//trio special handling: show UPD file is not empty
	if (type==GERMLINE_TRIO)
	{
		//UPDs
		FileLocation upd_loc = GlobalServiceProvider::fileLocationProvider().getAnalysisUpdFile();
		if (!upd_loc.exists)
		{
			QMessageBox::warning(this, "UPD detection", "The UPD file is missing!\n" + upd_loc.filename);
		}
		else
		{
			QStringList upd_data = Helper::loadTextFile(upd_loc.filename, false, QChar::Null, true);
			if (upd_data.count()>1)
			{
				ScrollableTextDialog dlg(this, "UPD(s) detected!");
				QStringList headers = upd_data[0].split("\t");
				for (int r=1; r<upd_data.count(); ++r)
				{
					QStringList parts = upd_data[r].split("\t");
					QString line = parts[0] + ":" + parts[1] + "-" + parts[2];
					for(int c=3 ; c<parts.count(); ++c)
					{
						line += " " + headers[c] + "=" + parts[c];
					}
					dlg.appendLine(line);
				}
				dlg.exec();
			}
		}
	}

	//check report sample ROH file exists
	QStringList roh_files = GlobalServiceProvider::fileLocationProvider().getRohFiles(false).filterById(germlineReportSample()).asStringList();
	if (roh_files.isEmpty())
	{
		QMessageBox::warning(this, "Runs of homozygosity", "Could not find a ROH file for sample " + germlineReportSample() + ". Aborting!");
		return;
	}

	RohWidget* list = new RohWidget(roh_files[0], ui_.filters);
	auto dlg = GUIHelper::createDialog(list, "Runs of homozygosity of " + variants_.analysisName());
	addModelessDialog(dlg);
}

void MainWindow::on_actionGeneSelector_triggered()
{
	if (filename_=="") return;
	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE) return;

	QString ps_name = germlineReportSample();

	//show dialog
	GeneSelectorDialog dlg(ps_name, this);
	if (dlg.exec())
	{
		//copy report to clipboard
		QString report = dlg.report();
		QApplication::clipboard()->setText(report);

		//show message
		if (QMessageBox::question(this, "Gene selection report", "Gene selection report was copied to clipboard.\nDo you want to open the sub-panel design dialog for selected genes?")==QMessageBox::Yes)
		{
			openSubpanelDesignDialog(dlg.genesForVariants());
		}
	}
}

void MainWindow::on_actionCircos_triggered()
{
	if (filename_=="") return;

	//load plot file
	FileLocationList plot_files = GlobalServiceProvider::fileLocationProvider().getCircosPlotFiles(true);
	if (plot_files.isEmpty()) return; //this should not happen because the button is not enabled then...
	if (!plot_files[0].exists)
	{
		QMessageBox::warning(this, "Circos plot file access", "Circos plot image file does not exist or the URL has expired");
		return;
	}
	//show plot
	CircosPlotWidget* widget = new CircosPlotWidget(plot_files[0].filename);
	auto dlg = GUIHelper::createDialog(widget, "Circos Plot of " + variants_.analysisName());
	addModelessDialog(dlg);
}

void MainWindow::on_actionExpressionData_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;

	QString title = "Expression data";

	NGSD db;
	QString sample_id = db.sampleId(filename_, false);
	if (sample_id=="")
	{
		QMessageBox::warning(this, title, "Error: Sample not found in NGSD!");
		return;
	}

	//get count files of all RNA processed samples corresponding to the current sample
	QStringList rna_ps_ids;
	foreach (int rna_sample, db.relatedSamples(sample_id.toInt(), "same sample", "RNA"))
	{
		rna_ps_ids << db.getValues("SELECT id FROM processed_sample WHERE sample_id=:0", QString::number(rna_sample));
	}

	QStringList rna_count_files;
	foreach (QString rna_ps_id, rna_ps_ids)
	{
		FileLocation file_location = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::EXPRESSION);
		if (file_location.exists) rna_count_files << file_location.filename;
	}
	rna_count_files.removeDuplicates();

	if (rna_count_files.isEmpty())
	{
		QMessageBox::warning(this, title, "Error: No RNA count files of corresponding RNA samples found!");
		return;
	}

	//select file to open
	QString count_file;
	if (rna_count_files.size()==1)
	{
		count_file = rna_count_files.at(0);
	}
	else
	{
		bool ok;
		count_file = getFileSelectionItem(title, "Multiple RNA count files found.\nPlease select a file:", rna_count_files, &ok);
		if (!ok) return;
	}

	int rna_sys_id = db.processingSystemIdFromProcessedSample(count_file);
	QString rna_ps_id = db.processedSampleId(count_file);
	QString tissue = db.getSampleData(db.sampleId(count_file)).tissue;
	QString project = db.getProcessedSampleData(rna_ps_id).project_name;

	GeneSet variant_target_region;
	if(ui_.filters->phenotypes().count() > 0)
	{
		foreach (const Phenotype& phenotype, ui_.filters->phenotypes())
		{
			variant_target_region << db.phenotypeToGenes(db.phenotypeIdByAccession(phenotype.accession()), false);
		}
	}

	if(ui_.filters->targetRegion().isValid())
	{
		if (variant_target_region.isEmpty())
		{
			variant_target_region = ui_.filters->targetRegion().genes;
		}
		else
		{
			variant_target_region = ui_.filters->targetRegion().genes.intersect(variant_target_region);
		}
	}

	RnaCohortDeterminationStategy cohort_type;
	if (germlineReportSupported())
	{
		cohort_type = RNA_COHORT_GERMLINE;
	}
	else
	{
		cohort_type = RNA_COHORT_SOMATIC;
	}

	ExpressionGeneWidget* widget = new ExpressionGeneWidget(count_file, rna_sys_id, tissue, ui_.filters->genes().toStringList().join(", "), variant_target_region, project, rna_ps_id,
															cohort_type, this);
	auto dlg = GUIHelper::createDialog(widget, "Gene expression of " + db.processedSampleName(rna_ps_id) + " (DNA: " + variants_.analysisName() + ")");
	addModelessDialog(dlg);
}

void MainWindow::on_actionExonExpressionData_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;

	QString title = "Exon expression data";

	NGSD db;
	QString sample_id = db.sampleId(filename_, false);
	if (sample_id=="")
	{
		QMessageBox::warning(this, title, "Error: Sample not found in NGSD!");
		return;
	}

	//get count files of all RNA processed samples corresponding to the current sample
	QStringList rna_ps_ids;
	foreach (int rna_sample, db.relatedSamples(sample_id.toInt(), "same sample", "RNA"))
	{
		rna_ps_ids << db.getValues("SELECT id FROM processed_sample WHERE sample_id=:0", QString::number(rna_sample));
	}

	QStringList rna_count_files;
	foreach (QString rna_ps_id, rna_ps_ids)
	{
		FileLocation file_location = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::EXPRESSION_EXON);
		if (file_location.exists) rna_count_files << file_location.filename;
	}
	rna_count_files.removeDuplicates();

	if (rna_count_files.isEmpty())
	{
		QMessageBox::warning(this, title, "Error: No RNA count files of corresponding RNA samples found!");
		return;
	}

	//select file to open
	QString count_file;
	if (rna_count_files.size()==1)
	{
		count_file = rna_count_files.at(0);
	}
	else
	{
		bool ok;
		count_file = getFileSelectionItem(title,"Multiple RNA count files found.\nPlease select a file:", rna_count_files, &ok);
		if (!ok) return;
	}

	int rna_sys_id = db.processingSystemIdFromProcessedSample(count_file);
	QString rna_ps_id = db.processedSampleId(count_file);
	QString tissue = db.getSampleData(db.sampleId(count_file)).tissue;
	QString project = db.getProcessedSampleData(rna_ps_id).project_name;

	GeneSet variant_target_region;
	if(ui_.filters->phenotypes().count() > 0)
	{
		foreach (const Phenotype& phenotype, ui_.filters->phenotypes())
		{
			variant_target_region << db.phenotypeToGenes(db.phenotypeIdByAccession(phenotype.accession()), false);
		}
	}

	if(ui_.filters->targetRegion().isValid())
	{
		if (variant_target_region.isEmpty())
		{
			variant_target_region = ui_.filters->targetRegion().genes;
		}
		else
		{
			variant_target_region = ui_.filters->targetRegion().genes.intersect(variant_target_region);
		}
	}

	RnaCohortDeterminationStategy cohort_type = RNA_COHORT_GERMLINE;
	if (somaticReportSupported()) cohort_type = RNA_COHORT_SOMATIC;

	ExpressionExonWidget* widget = new ExpressionExonWidget(count_file, rna_sys_id, tissue, ui_.filters->genes().toStringList().join(", "), variant_target_region, project, rna_ps_id, cohort_type, this);
	auto dlg = GUIHelper::createDialog(widget, "Exon expression of " + db.processedSampleName(rna_ps_id) + " (DNA: " + variants_.analysisName() + ")");
	addModelessDialog(dlg);
}

void MainWindow::on_actionShowSplicing_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;

	NGSD db;

	//get all available files
	QStringList splicing_files;
	foreach (int rna_sample_id, db.relatedSamples(db.sampleId(variants_.mainSampleName()).toInt(), "same sample", "RNA"))
	{
		// check for required files
		foreach (const QString& rna_ps_id, db.getValues("SELECT id FROM processed_sample WHERE sample_id=:0", QString::number(rna_sample_id)))
		{
			// search for fusion file
			FileLocation splicing_file = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::SPLICING_ANN);
			if (splicing_file.exists) splicing_files << splicing_file.filename;
		}
	}

	if (splicing_files.isEmpty())
	{
		QMessageBox::warning(this, "Splicing files missing", "Error: No RNA splicing files of corresponding RNA samples found!");
		return;
	}

	//select file to open
	QString splicing_filepath;
	if (splicing_files.size()==1)
	{
		splicing_filepath = splicing_files.at(0);
	}
	else
	{
		bool ok;
		splicing_filepath = getFileSelectionItem("Multiple files found", "Multiple RNA splicing files found.\nPlease select a file:", splicing_files, &ok);
		if (!ok) return;
	}

	SplicingWidget* splicing_widget = new SplicingWidget(splicing_filepath, this);

	auto dlg = GUIHelper::createDialog(splicing_widget, "Splicing Alterations of " + variants_.analysisName());
	addModelessDialog(dlg);
}

void MainWindow::on_actionShowRnaFusions_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;

	NGSD db;

	//get all available files
	QStringList arriba_fusion_files;
	foreach (int rna_sample_id, db.relatedSamples(db.sampleId(variants_.mainSampleName()).toInt(), "same sample", "RNA"))
	{
		// check for required files
		foreach (const QString& rna_ps_id, db.getValues("SELECT id FROM processed_sample WHERE sample_id=:0", QString::number(rna_sample_id)))
		{
			// search for fusion file
			FileLocation fusion_file = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::FUSIONS);
			if (fusion_file.exists) arriba_fusion_files << fusion_file.filename;
		}
	}

	if (arriba_fusion_files.isEmpty())
	{
		QMessageBox::warning(this, "Fusion files missing", "Error: No RNA fusion files of corresponding RNA samples found!");
		return;
	}

	//select file to open
	QString fusion_filepath;
	if (arriba_fusion_files.size()==1)
	{
		fusion_filepath = arriba_fusion_files.at(0);
	}
	else
	{
		bool ok;
		fusion_filepath = getFileSelectionItem("Multiple files found", "Multiple RNA fusion files found.\nPlease select a file:", arriba_fusion_files, &ok);
		if (!ok) return;
	}

	FusionWidget* fusion_widget = new FusionWidget(fusion_filepath, this);

	auto dlg = GUIHelper::createDialog(fusion_widget, "Fusions of " + variants_.analysisName() + " (arriba)");
	addModelessDialog(dlg);
}

void MainWindow::on_actionShowProcessingSystemCoverage_triggered()
{
	//set filter widget
	FilterWidget* variant_filter_widget = nullptr;
	if(filename_ != "") variant_filter_widget = ui_.filters;

	auto expression_level_widget = new ExpressionOverviewWidget(variant_filter_widget, this);

	auto dlg = GUIHelper::createDialog(expression_level_widget, "Expression of processing systems");
	addModelessDialog(dlg);
}

void MainWindow::on_actionRE_triggered()
{
	if (filename_=="") return;
	if (variants_.type()!=GERMLINE_SINGLESAMPLE) return;

	// determine repeat expansion file name
	FileLocationList re_files = GlobalServiceProvider::fileLocationProvider().getRepeatExpansionFiles(false);
	if (re_files.isEmpty()) return; //this should not happen because the button is not enabled then...

	QString ps_name = variants_.mainSampleName();

	//get sample type
	bool is_exome = false;
	if (LoginManager::active())
	{
		NGSD db;
		QString ps_id = db.processedSampleId(ps_name, false);
		is_exome = ps_id!="" && db.getProcessedSampleData(ps_id).processing_system_type=="WES";
	}

	//show dialog
	RepeatExpansionWidget* widget = new RepeatExpansionWidget(re_files[0].filename, is_exome);
	auto dlg = GUIHelper::createDialog(widget, "Repeat Expansions of " + variants_.analysisName());

	addModelessDialog(dlg);
}

void MainWindow::on_actionPRS_triggered()
{
	if (filename_=="") return;

	// determine PRS file name
	FileLocationList prs_files = GlobalServiceProvider::fileLocationProvider().getPrsFiles(false);
	if (prs_files.isEmpty()) return; //this should not happen because the button is not enabled then...

	//show dialog
	PRSWidget* widget = new PRSWidget(prs_files[0].filename);
	auto dlg = GUIHelper::createDialog(widget, "Polygenic Risk Scores of " + variants_.analysisName());
	addModelessDialog(dlg);
}

void MainWindow::on_actionDesignCfDNAPanel_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;
	if (!(somaticReportSupported()||tumoronlyReportSupported())) return;

	// Workaround to manual add panels for non patient-specific processing systems
	DBTable cfdna_processing_systems = NGSD().createTable("processing_system", "SELECT id, name_short FROM processing_system WHERE type='cfDNA (patient-specific)' OR type='cfDNA'");
	// TODO: reactivate
//	DBTable cfdna_processing_systems = NGSD().createTable("processing_system", "SELECT id, name_short FROM processing_system WHERE type='cfDNA (patient-specific)'");

	QSharedPointer<CfDNAPanelDesignDialog> dialog(new CfDNAPanelDesignDialog(variants_, filter_result_, somatic_report_settings_.report_config, variants_.mainSampleName(), cfdna_processing_systems, this));
	dialog->setWindowFlags(Qt::Window);

	addModelessDialog(dialog);
}

void MainWindow::on_actionShowCfDNAPanel_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;
	if (!(somaticReportSupported()||tumoronlyReportSupported())) return;

	NGSD db;
// get cfDNA panels:
	QList<CfdnaPanelInfo> cfdna_panels = db.cfdnaPanelInfo(db.processedSampleId(variants_.mainSampleName()));
	CfdnaPanelInfo selected_panel;
	if (cfdna_panels.size() < 1)	{
		// show message
		GUIHelper::showMessage("No cfDNA panel found!", "No cfDNA panel was found for the given tumor sample!");
		return;
	}
	else if (cfdna_panels.size() > 1)
	{
		QStringList cfdna_panel_description;
		foreach (const CfdnaPanelInfo& panel, cfdna_panels)
		{
			cfdna_panel_description.append("cfDNA panel for " + db.getProcessingSystemData(panel.processing_system_id).name  + " (" + panel.created_date.toString("dd.MM.yyyy") + " by "
										   + db.userName(panel.created_by) + ")");
		}

		QComboBox* cfdna_panel_selector = new QComboBox(this);
		cfdna_panel_selector->addItems(cfdna_panel_description);

		// create dlg
		auto dlg = GUIHelper::createDialog(cfdna_panel_selector, "Select cfDNA panel", "", true);
		int btn = dlg->exec();
		if (btn!=1)
		{
			return;
		}
		selected_panel = cfdna_panels.at(cfdna_panel_selector->currentIndex());
	}
	else
	{
		selected_panel = cfdna_panels.at(0);
	}

	//show dialog
	CfDNAPanelWidget* widget = new CfDNAPanelWidget(selected_panel);
	auto dlg = GUIHelper::createDialog(widget, "cfDNA panel for tumor " + variants_.analysisName());
	addModelessDialog(dlg);
}

void MainWindow::on_actionCfDNADiseaseCourse_triggered()
{
	if (filename_=="") return;
	if (!somaticReportSupported()) return;

	DiseaseCourseWidget* widget = new DiseaseCourseWidget(variants_.mainSampleName());
	auto dlg = GUIHelper::createDialog(widget, "Personalized cfDNA variants");
	addModelessDialog(dlg);
}

void MainWindow::on_actionCfDNAAddExcludedRegions_triggered()
{
	if (filename_=="") return;
	if (!LoginManager::active()) return;
	if (!somaticReportSupported()) return;

	QSharedPointer<cfDNARemovedRegions> dialog(new cfDNARemovedRegions(variants_.mainSampleName(), this));
	dialog->setWindowFlags(Qt::Window);

	addModelessDialog(dialog);
}

void MainWindow::on_actionGeneOmimInfo_triggered()
{
	GeneOmimInfoWidget* widget = new GeneOmimInfoWidget(this);
	auto dlg = GUIHelper::createDialog(widget, "OMIM information for genes");
	dlg->exec();
}

void MainWindow::openVariantListFolder()
{
	if (filename_=="") return;

	try
	{
		QString sample_folder;

		if (GlobalServiceProvider::fileLocationProvider().isLocal())
		{
			sample_folder = QFileInfo(filename_).absolutePath();
		}
		else //client-server (allow opening folders if GSvar project paths are configured)
		{
			NGSD db;
			if (variants_.type()==AnalysisType::GERMLINE_SINGLESAMPLE)
			{
				QString ps = variants_.mainSampleName();
				QString ps_id = db.processedSampleId(ps);
				QString project_type = db.getProcessedSampleData(ps_id).project_type;
				QString project_folder = db.projectFolder(project_type).trimmed();
				if (!project_folder.isEmpty())
				{
					sample_folder = db.processedSamplePath(ps_id, PathType::SAMPLE_FOLDER);
					if (!QDir(sample_folder).exists()) THROW(Exception, "Sample folder does not exist: " + sample_folder);
				}
			}
			else
			{
				THROW(Exception, "In client-server mode, opening analysis folders is only supported for germline single sample!");
			}
		}

		QDesktopServices::openUrl(sample_folder);
	}
	catch(Exception& e)
	{
		QMessageBox::information(this, "Open analysis folder", "Could not open analysis folder:\n" + e.message());
	}
}

void MainWindow::openVariantListQcFiles()
{
	if (filename_=="") return;

	const FileLocationProvider& flp = GlobalServiceProvider::fileLocationProvider();

	foreach(const FileLocation& file, flp.getQcFiles())
	{
		if (flp.isLocal())
		{
			QDesktopServices::openUrl(QUrl::fromLocalFile(file.filename));
		}
		else
		{
			QDesktopServices::openUrl(file.filename);
		}
	}
}

void MainWindow::on_actionReanalyze_triggered()
{
	if (filename_=="") return;

	AnalysisType type = variants_.type();
	SampleHeaderInfo header_info = variants_.getSampleHeader();
	QList<AnalysisJobSample> samples;
	if (type==GERMLINE_SINGLESAMPLE  || type==CFDNA || type==SOMATIC_SINGLESAMPLE)
	{
		samples << AnalysisJobSample {header_info[0].id, ""};
	}
	else if (type==GERMLINE_MULTISAMPLE)
	{
		foreach(const SampleInfo& info, header_info)
		{
			samples << AnalysisJobSample {info.id, info.isAffected() ? "affected" : "control"};
		}
	}
	else if (type==GERMLINE_TRIO)
	{
		foreach(const SampleInfo& info, header_info)
		{
			if(info.isAffected())
			{
				samples << AnalysisJobSample {info.id, "child"};
			}
			else
			{
				samples << AnalysisJobSample {info.id, info.gender()=="male" ? "father" : "mother"};
			}
		}
	}
	else if (type==SOMATIC_PAIR)
	{
		foreach(const SampleInfo& info, header_info)
		{
			samples << AnalysisJobSample {info.id, info.isTumor() ? "tumor" : "normal"};
		}
	}

	GSvarHelper::queueSampleAnalysis(type, samples, this);
}

void MainWindow::delayedInitialization()
{
	//initialize LOG file
	Log::setFileEnabled(true);
	Log::appInfo();

	//check that INI file is configured
	QStringList keys_missing;
	foreach(QString key, QStringList() << "build" << "reference_genome" << "igv_app" << "igv_genome" << "threads")
	{
	   if (!Settings::contains(key)) keys_missing << key;
	}
	if (!keys_missing.isEmpty())
	{
		QMessageBox::warning(this, "GSvar setup error", "The GSvar INI file is not set up correctly.\nThe following keys are missing or contain no value: " + keys_missing.join(", ") + "\nPlease inform your administrator!");
		close();
		return;
	}

	//close the app if the server is not available
	if (ClientHelper::isClientServerMode())
	{
		if (!isServerRunning())
		{
			close();
			return;
		}
	}

	//user login for database
	if (GlobalServiceProvider::database().enabled())
	{
		LoginDialog dlg(this);
		dlg.exec();

		if (LoginManager::active())
		{
			try
			{
				ui_.filters->loadTargetRegions();
			}
			catch(Exception& e)
			{
				Log::warn("Target region data for filter widget could not be loaded from NGSD: " + e.message());
			}
		}
	}

	//create default IGV session (variants)
	IGVSession* igv_default = IgvSessionManager::create(this, "Default IGV", Settings::path("igv_app").trimmed(), Settings::string("igv_host"), Settings::path("igv_genome"));
	connect(igv_default, SIGNAL(started()), this, SLOT(changeIgvIconToActive()));
	connect(igv_default, SIGNAL(finished()), this, SLOT(changeIgvIconToNormal()));

	//create IGV session for virus detection
	QString virus_genome;
	if (!ClientHelper::isClientServerMode())
	{
		virus_genome = Settings::string("igv_virus_genome", true);
		if (virus_genome.isEmpty()) QMessageBox::information(this, "Virus genome not set", "Virus genome path is missing from the settings!");
	}
	else
	{
	   virus_genome = ClientHelper::serverApiUrl() + "genome/somatic_viral.fa";
	}
	IGVSession* igv_virus = IgvSessionManager::create(this, "Virus IGV", Settings::path("igv_app").trimmed(), Settings::string("igv_host"), virus_genome);
	connect(igv_virus, SIGNAL(started()), this, SLOT(changeIgvIconToActive()));
	connect(igv_virus, SIGNAL(finished()), this, SLOT(changeIgvIconToNormal()));

	//init GUI
	updateRecentSampleMenu();
	updateIGVMenu();
	updateNGSDSupport();
	registerCustomContextMenuActions();

	//parse arguments
	for (int i=1; i<QApplication::arguments().count(); ++i)
	{
		QString arg = QApplication::arguments().at(i);

		if (i==1) //first argument: sample
		{
			if (QFile::exists(arg)) //file path
			{
				loadFile(arg);
			}
			else if (LoginManager::active()) //processed sample name (via NGSD)
			{
				NGSD db;
				if (db.processedSampleId(arg, false)!="")
				{
					openProcessedSampleFromNGSD(arg, false);
				}
				else if (db.sampleId(arg, false)!="")
				{
					openSampleFromNGSD(arg);
				}
			}
		}
		else if (arg.startsWith("filter:")) //filter (by name)
		{
			int sep_pos = arg.indexOf(':');
			QString filter_name = arg.mid(sep_pos+1).trimmed();

			if (!ui_.filters->setFilter(filter_name))
			{
				qDebug() << "Filter name " << arg << " not found. Ignoring it!";
			}
		}
		else if (arg.startsWith("roi:")) //target region (by name)
		{
			int sep_pos = arg.indexOf(':');
			QString roi_name = arg.mid(sep_pos+1).trimmed();

			if (!ui_.filters->setTargetRegionByDisplayName(roi_name))
			{
				qDebug() << "Target region name " << roi_name << " not found. Ignoring it!";
			}
		}
		else
		{
			qDebug() << "Unprocessed argument: " << arg;
		}
	}
}

void MainWindow::variantCellDoubleClicked(int row, int /*col*/)
{
    const Variant& v = variants_[ui_.vars->rowToVariantIndex(row)];
    IgvSessionManager::get(0).gotoInIGV(v.chr().str() + ":" + QString::number(v.start()) + "-" + QString::number(v.end()), true);
}

void MainWindow::variantHeaderDoubleClicked(int row)
{
	if (!LoginManager::active()) return;

	int var_index = ui_.vars->rowToVariantIndex(row);
	editVariantReportConfiguration(var_index);
}

void MainWindow::openCustomIgvTrack()
{
	QAction* action = qobject_cast<QAction*>(sender());
	if (action==nullptr) return;

    IgvSessionManager::get(0).loadFileInIGV(action->toolTip(), false);
}

void MainWindow::editVariantValidation(int index)
{
	Variant& variant = variants_[index];

	try
	{
		QString ps = selectProcessedSample();
		if (ps.isEmpty()) return;

		NGSD db;

		//get variant ID - add if missing
		QString variant_id = db.variantId(variant, false);
		if (variant_id=="")
		{
			variant_id = db.addVariant(variant, variants_);
		}

		//get sample ID
		QString sample_id = db.sampleId(ps);

		//get variant validation ID - add if missing
		QVariant val_id = db.getValue("SELECT id FROM variant_validation WHERE variant_id='" + variant_id + "' AND sample_id='" + sample_id + "'", true);
		bool added_validation_entry = false;
		if (!val_id.isValid())
		{
			//get genotype
			QByteArray genotype = "het";
			int i_genotype = variants_.getSampleHeader().infoByID(ps).column_index;
			if (i_genotype!=-1) //genotype column only available in germline, but not for somatic analysis.
			{
				genotype = variant.annotations()[i_genotype];
			}

			//insert
			SqlQuery query = db.getQuery();
			query.exec("INSERT INTO variant_validation (user_id, sample_id, variant_type, variant_id, genotype, status) VALUES ('" + LoginManager::userIdAsString() + "','" + sample_id + "','SNV_INDEL','" + variant_id + "','" + genotype + "','n/a')");
			val_id = query.lastInsertId();

			added_validation_entry = true;
		}

		ValidationDialog dlg(this, val_id.toInt());

		if (dlg.exec())
		{
			//update DB
			dlg.store();

			//update variant table
			QByteArray status = dlg.status().toUtf8();
			if (status=="true positive") status = "TP";
			if (status=="false positive") status = "FP";
			int i_validation = variants_.annotationIndexByName("validation", true, true);
			variant.annotations()[i_validation] = status;

			//update details widget and filtering
			ui_.variant_details->updateVariant(variants_, index);
			refreshVariantTable();

			//mark variant list as changed
			markVariantListChanged(variant, "validation", status);
		}
		else if (added_validation_entry)
		{
			// remove created but empty validation if ValidationDialog is aborted
			SqlQuery query = db.getQuery();
			query.exec("DELETE FROM variant_validation WHERE id=" + val_id.toString());
		}
	}
	catch (DatabaseException& e)
	{
		GUIHelper::showMessage("NGSD error", e.message());
	}
}

void MainWindow::editVariantComment(int index)
{
	Variant& variant = variants_[index];

	try
	{
		//add variant if missing
		NGSD db;
		if (db.variantId(variant, false)=="")
		{
			db.addVariant(variant, variants_);
		}

		bool ok = true;
		QByteArray text = QInputDialog::getMultiLineText(this, "Variant comment", "Text: ", db.comment(variant), &ok).toUtf8();

		if (ok)
		{
			//update DB
			db.setComment(variant, text);

			//update datastructure (if comment column is present)
			int col_index = variants_.annotationIndexByName("comment", true, false);
			if (col_index!=-1)
			{
				variant.annotations()[col_index] = text;
				refreshVariantTable();

				//mark variant list as changed
				markVariantListChanged(variant, "comment", text);
			}
		}
	}
	catch (DatabaseException& e)
	{
		GUIHelper::showMessage("NGSD error", e.message());
	}
}

void MainWindow::showAfHistogram_all()
{
	showAfHistogram(false);
}

void MainWindow::showAfHistogram_filtered()
{
	showAfHistogram(true);
}

void MainWindow::showCnHistogram()
{
	if (filename_=="") return;

	QString title = "Copy-number histogram";

	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE)
	{
		QMessageBox::information(this, title, "This functionality is only available for germline single sample and germline trio analysis.");
		return;
	}

	//check report sample SEG file exists
	QStringList seg_files = GlobalServiceProvider::fileLocationProvider().getCnvCoverageFiles(false).filterById(germlineReportSample()).asStringList();
	if (seg_files.isEmpty())
	{
		QMessageBox::warning(this, title, "Could not find a SEG file for sample " + germlineReportSample() + ". Aborting!");
		return;
	}

	try
	{
		//get region
		Chromosome chr;
		int start, end;
		QString region_text = QInputDialog::getText(this, title, "genomic region");
		if (region_text=="") return;

		NGSHelper::parseRegion(region_text, chr, start, end, true);

		//determine CN values
		QVector<double> cn_values;
		VersatileTextStream stream(seg_files[0]);
		while (!stream.atEnd())
		{
			QString line = stream.readLine();
			QStringList parts = line.split("\t");
			if (parts.count()<6) continue;

			//check if range overlaps input interval
			Chromosome chr2(parts[1]);
			if (chr!=chr2) continue;

			int start2 = Helper::toInt(parts[2], "Start coordinate");
			int end2 = Helper::toInt(parts[3], "End coordinate");
			if (!BasicStatistics::rangeOverlaps(start, end, start2, end2)) continue;

			//skip invalid copy-numbers
			QString cn_str = parts[5];
			if (cn_str.toLower()=="nan") continue;
			double cn = Helper::toDouble(cn_str, "Copy-number");
			if (cn<0) continue;

			cn_values << cn;
		}

		//create histogram
		std::sort(cn_values.begin(), cn_values.end());
		double median = BasicStatistics::median(cn_values,false);
		double max = ceil(median*2+0.0001);
		Histogram hist(0.0, max, max/40);
		foreach(double cn, cn_values)
		{
			hist.inc(cn, true);
		}

		//show chart
		QChartView* view = GUIHelper::histogramChart(hist, "Copy-number");
		auto dlg = GUIHelper::createDialog(view, QString("Copy-number in region ") + region_text);
		dlg->exec();
	}
	catch(Exception& e)
	{
		QMessageBox::warning(this, title, "Error:\n" + e.message());
		return;
	}
}

void MainWindow::showBafHistogram()
{
	if (filename_=="") return;

	QString title = "BAF histogram";

	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE)
	{
		QMessageBox::information(this, title, "This functionality is only available for germline single sample and germline trio analysis.");
		return;
	}

	//check report sample SEG file exists
	QStringList baf_files = GlobalServiceProvider::fileLocationProvider().getBafFiles(false).filterById(germlineReportSample()).asStringList();
	if (baf_files.isEmpty())
	{
		QMessageBox::warning(this, title, "Could not find a BAF file for sample " + germlineReportSample() + ". Aborting!");
		return;
	}

	try
	{
		//get region
		Chromosome chr;
		int start, end;
		QString region_text = QInputDialog::getText(this, title, "genomic region");
		if (region_text=="") return;

		NGSHelper::parseRegion(region_text, chr, start, end, true);

		//determine CN values
		Histogram hist(0.0, 1.0, 0.025);
        QSharedPointer<VersatileFile> file = Helper::openVersatileFileForReading(baf_files[0]);
        while (!file->atEnd())
        {
            QString line = file->readLine();
			QStringList parts = line.split("\t");
			if (parts.count()<5) continue;

			//check if range overlaps input interval
			Chromosome chr2(parts[0]);
			if (chr!=chr2) continue;

			int start2 = Helper::toInt(parts[1], "Start coordinate");
			int end2 = Helper::toInt(parts[2], "End coordinate");
			if (!BasicStatistics::rangeOverlaps(start, end, start2, end2)) continue;

			double baf =  Helper::toDouble(parts[4], "BAF");
			hist.inc(baf, true);
		}

		//show chart
		QChartView* view = GUIHelper::histogramChart(hist, "BAF");
		auto dlg = GUIHelper::createDialog(view, QString("BAF in region ") + region_text);
		dlg->exec();
	}
	catch(Exception& e)
	{
		QMessageBox::warning(this, title, "Error:\n" + e.message());
		return;
	}
}

void MainWindow::showAfHistogram(bool filtered)
{
	if (filename_=="") return;

	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO)
	{
		QMessageBox::information(this, "Allele frequency histogram", "This functionality is only available for germline single sample and germline trio analysis.");
		return;
	}

	//create histogram
	Histogram hist(0.0, 1.0, 0.05);
	int col_quality = variants_.annotationIndexByName("quality");
	for (int i=0; i<variants_.count(); ++i)
	{

		if (filtered && !filter_result_.passing(i)) continue;

		QByteArrayList parts = variants_[i].annotations()[col_quality].split(';');
		foreach(const QByteArray& part, parts)
		{
			if (part.startsWith("AF="))
			{
				bool ok;
				QString value_str = part.mid(3);
				if (type==GERMLINE_TRIO) value_str = value_str.split(',')[0];
				double value = value_str.toDouble(&ok);
				if (ok)
				{
					hist.inc(value, true);
				}
			}
		}
	}

	//show chart
	QChartView* view = GUIHelper::histogramChart(hist, "Allele frequency");
	auto dlg = GUIHelper::createDialog(view, QString("Allele frequency histogram ") + (filtered ? "(after filter)" : "(all variants)"));
	dlg->exec();
}

void MainWindow::on_actionEncrypt_triggered()
{
	//get input
	QString input = QInputDialog::getText(this, "Text for encryption", "text");
	if (input.isEmpty()) return;

	//decrypt
	QStringList out_lines;
	out_lines << ("Input text: " + input);
	try
	{
		qulonglong crypt_key = ToolBase::encryptionKey("encryption helper");
		out_lines << ("Encrypted text: " + SimpleCrypt(crypt_key).encryptToString(input));
	}
	catch(Exception& e)
	{
		out_lines << ("Error: " + e.message());
	}

	//show output
	QTextEdit* edit = new QTextEdit(this);
	edit->setText(out_lines.join("\n"));
	auto dlg = GUIHelper::createDialog(edit, "Encryption output");
	dlg->exec();
}

void MainWindow::on_actionSettings_triggered()
{
	SettingsDialog dlg(this);
	if (dlg.exec()==QDialog::Accepted)
	{
		dlg.storeSettings();
	}
}

void MainWindow::on_actionSampleSearch_triggered()
{
	SampleSearchWidget* widget = new SampleSearchWidget(this);
	openTab(QIcon(":/Icons/NGSD_sample_search.png"), "Sample search", widget);
}

void MainWindow::on_actionRunOverview_triggered()
{
	SequencingRunOverview* widget = new SequencingRunOverview(this);
	openTab(QIcon(":/Icons/NGSD_run_overview.png"), "Sequencing run overview", widget);
}

void MainWindow::addModelessDialog(QSharedPointer<QDialog> dlg, bool maximize)
{
	if (maximize)
	{
		dlg->showMaximized();
	}
	else
	{
		dlg->show();
	}

	connect(dlg.data(), SIGNAL(finished(int)), this, SLOT(deleteClosedModelessDialogs()));

	modeless_dialogs_.append(dlg);
}

void MainWindow::deleteClosedModelessDialogs()
{
	//Clean up when we add another dialog. Like that, only one dialog can be closed and not destroyed at the same time.
	for (int i=modeless_dialogs_.count()-1; i>=0; --i)
	{
		if (modeless_dialogs_[i]->isHidden())
		{
			modeless_dialogs_.removeAt(i);
		}
	}
}

void MainWindow::importPhenotypesFromNGSD()
{
	QString ps_name = germlineReportSupported() ? germlineReportSample() : variants_.mainSampleName();
	try
	{
		NGSD db;
		QString sample_id = db.sampleId(ps_name);
		PhenotypeList phenotypes = db.getSampleData(sample_id).phenotypes;

		ui_.filters->setPhenotypes(phenotypes);
	}
	catch(Exception& e)
	{
		GUIHelper::showException(this, e, "Error loading phenotype data from NGSD");
	}
}

void MainWindow::createSubPanelFromPhenotypeFilter()
{
	//convert phenotypes to genes
	QApplication::setOverrideCursor(Qt::BusyCursor);
	NGSD db;
	GeneSet genes;
	foreach(const Phenotype& pheno, ui_.filters->phenotypes())
	{
		genes << db.phenotypeToGenes(db.phenotypeIdByAccession(pheno.accession()), true);
	}
	QApplication::restoreOverrideCursor();

	//open dialog
	openSubpanelDesignDialog(genes);
}

void MainWindow::on_actionOpen_triggered()
{
	//get file to open
	QString path = Settings::path("path_variantlists", true);
	QString filename = QFileDialog::getOpenFileName(this, "Open variant list", path, "GSvar files (*.GSvar);;All files (*.*)");
	if (filename=="") return;

	//update data
	loadFile(filename);
}

void MainWindow::on_actionOpenByName_triggered()
{
	ProcessedSampleSelector dlg(this, false);
	dlg.showSearchMulti(true);
	if (!dlg.exec()) return;

	QString ps_name = dlg.processedSampleName();
	if (ps_name.isEmpty()) return;
	openProcessedSampleFromNGSD(ps_name, dlg.searchMulti());
}

void MainWindow::openProcessedSampleFromNGSD(QString processed_sample_name, bool search_multi)
{
	checkClientUpdates();
	try
	{
		NGSD db;
		QString processed_sample_id = db.processedSampleId(processed_sample_name);

		//check user can access
		if (!db.userCanAccess(LoginManager::userId(), processed_sample_id.toInt()))
		{
			INFO(AccessDeniedException, "You do not have permissions to open sample '" + processed_sample_name + "'!");
		}

		//processed sample exists > add to recent samples menu
		addToRecentSamples(processed_sample_name);

		//germline single sample analysis
		QStringList analyses;
		FileLocation file_location = GlobalServiceProvider::database().processedSamplePath(processed_sample_id, PathType::GSVAR);
		if (file_location.exists) analyses << file_location.filename;

		//somatic tumor sample > ask user if he wants to open the tumor-normal pair
		QString normal_sample = db.normalSample(processed_sample_id);
		if (normal_sample!="")
		{
			analyses << GlobalServiceProvider::database().secondaryAnalyses(processed_sample_name + "-" + normal_sample, "somatic");
		}
		//check for germline trio/multi analyses
		else if (search_multi)
		{
			analyses << GlobalServiceProvider::database().secondaryAnalyses(processed_sample_name, "trio");
			analyses << GlobalServiceProvider::database().secondaryAnalyses(processed_sample_name, "multi sample");
		}

		//determine analysis to load
		QString file;
		if (analyses.count()==0)
		{
			INFO(ArgumentException, "The GSvar file does not exist:\n" + file_location.filename);
		}
		else if (analyses.count()==1)
		{
			file = analyses[0];
		}
		else //several analyses > let the user decide
		{
			//create list of anaylsis names
			QList<MultiSampleAnalysisInfo> analysis_info_list = GlobalServiceProvider::database().getMultiSampleAnalysisInfo(analyses);
			QStringList names;
			foreach(MultiSampleAnalysisInfo info, analysis_info_list)
			{
				names.append(info.analysis_name);
			}

			//show selection dialog (analysis name instead of file name)
			bool ok = false;
			QString name = QInputDialog::getItem(this, "Several analyses of the sample present", "select analysis:", names, 0, false, &ok);
			if (!ok) return;

			int index = names.indexOf(name);
			foreach(QString ps_id, analysis_info_list[index].ps_sample_id_list)
			{
				if (!db.userCanAccess(LoginManager::userId(), ps_id.toInt()))
				{
					INFO(AccessDeniedException, "You do not have permissions to open all the samples from the selected multi-sample analysis!");
				}
			}
			file = analysis_info_list[index].analysis_file;
		}

		loadFile(file);
	}
	catch (Exception& e)
	{
		GUIHelper::showException(this, e, "Error opening processed sample by name");
	}
}

void MainWindow::openSampleFromNGSD(QString sample_name)
{
	try
	{
		NGSD db;
		QStringList processed_samples = db.getValues("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND ps.quality!='bad' AND ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples) AND s.name=:0", sample_name);
		if (processed_samples.isEmpty())
		{
			THROW(ArgumentException, "No high-quality processed sample found for sample name '" + sample_name + "'");
		}

		if (processed_samples.count()==1)
		{
			openProcessedSampleFromNGSD(processed_samples[0], false);
		}
		else
		{
			bool ok = false;
			QString ps = QInputDialog::getItem(this, "Several processed samples found for sample '" + sample_name + "'", "select processed sample:", processed_samples, 0, false, &ok);
			if (!ok) return;

			openProcessedSampleFromNGSD(ps, false);
		}
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Error opening sample from NGSD", e.message());
	}
}

void MainWindow::checkMendelianErrorRate(double cutoff_perc)
{
	QString output = "";
	try
	{
		SampleHeaderInfo infos = variants_.getSampleHeader();
		int col_c = infos.infoByStatus(true).column_index;

		bool above_cutoff = false;
		QStringList mers;
		foreach(const SampleInfo& info, infos)
		{
			if (info.isAffected()) continue;

			int errors = 0;
			int autosomal = 0;

			int col_p = info.column_index;

			for (int i=0; i<variants_.count(); ++i)
			{
				const Variant& v = variants_[i];
				if (!v.chr().isAutosome()) continue;
				++autosomal;

				QString geno_c = v.annotations()[col_c];
				QString geno_p = v.annotations()[col_p];

				if ((geno_p=="hom" && geno_c=="wt") || (geno_p=="wt" && geno_c=="hom")) ++errors;
			}

			double percentage = 100.0 * errors / autosomal;
			if (percentage>cutoff_perc) above_cutoff = true;
			mers << infos.infoByStatus(true).id + " - " + info.id + ": " + QString::number(errors) + "/" + QString::number(autosomal) + " ~ " + QString::number(percentage, 'f', 2) + "%";
		}

		if (above_cutoff)
		{
			output = "Mendelian error rate too high:\n" + mers.join("\n");
		}
	}
	catch (Exception& e)
	{
		output = "Mendelian error rate calulation not possible:\n" + e.message();
	}

	if (!output.isEmpty())
	{
		QMessageBox::warning(this, "Medelian error rate", output);
	}
}

void MainWindow::openProcessedSampleTab(QString ps_name)
{
	try
	{
		QString ps_id = NGSD().processedSampleId(ps_name);

		ProcessedSampleWidget* widget = new ProcessedSampleWidget(this, ps_id);
		connect(widget, SIGNAL(clearMainTableSomReport(QString)), this, SLOT(clearSomaticReportSettings(QString)));
		connect(widget, SIGNAL(addModelessDialog(QSharedPointer<QDialog>, bool)), this, SLOT(addModelessDialog(QSharedPointer<QDialog>, bool)));
		int index = openTab(QIcon(":/Icons/NGSD_sample.png"), ps_name, widget);
		if (Settings::boolean("debug_mode_enabled"))
		{
			ui_.tabs->setTabToolTip(index, "NGSD ID: " + ps_id);
		}
	}
	catch (Exception& e)
	{
		GUIHelper::showException(this, e, "Open processed sample tab");
	}
}

void MainWindow::openRunTab(QString run_name)
{
	QString run_id;
	try
	{
		run_id = NGSD().getValue("SELECT id FROM sequencing_run WHERE name=:0", true, run_name).toString();
	}
	catch (DatabaseException e)
	{
		GUIHelper::showMessage("NGSD error", "The run database ID could not be determined for '"  + run_name + "'!\nError message: " + e.message());
		return;
	}

	SequencingRunWidget* widget = new SequencingRunWidget(this, run_id);
	connect(widget, SIGNAL(addModelessDialog(QSharedPointer<QDialog>, bool)), this, SLOT(addModelessDialog(QSharedPointer<QDialog>, bool)));
	int index = openTab(QIcon(":/Icons/NGSD_run.png"), run_name, widget);
	if (Settings::boolean("debug_mode_enabled"))
	{
		ui_.tabs->setTabToolTip(index, "NGSD ID: " + run_id);
	}
}

void MainWindow::openGeneTab(QString symbol)
{
	QPair<QString, QString> approved = NGSD().geneToApprovedWithMessage(symbol);
	if (approved.second.startsWith("ERROR:"))
	{
		GUIHelper::showMessage("NGSD error", "Gene name '" + symbol + "' is not a HGNC-approved name!\nError message:\n" + approved.second);
		return;
	}

	GeneWidget* widget = new GeneWidget(this, symbol.toUtf8());
	int index = openTab(QIcon(":/Icons/NGSD_gene.png"), symbol, widget);
	if (Settings::boolean("debug_mode_enabled"))
	{
		ui_.tabs->setTabToolTip(index, "NGSD ID: " + QString::number(NGSD().geneId(symbol.toUtf8())));
	}
}

void MainWindow::openVariantTab(Variant variant)
{
	try
	{
		//check variant is in NGSD
		NGSD db;
		QString v_id = db.variantId(variant);

		//open tab
		VariantWidget* widget = new VariantWidget(variant, this);
		int index = openTab(QIcon(":/Icons/NGSD_variant.png"), variant.toString(), widget);

		//add database id
		if (Settings::boolean("debug_mode_enabled"))
		{
			ui_.tabs->setTabToolTip(index, "NGSD ID: " + v_id);
		}
	}
	catch(Exception& e)
	{
		QMessageBox::warning(this, "Open variant tab", e.message());
	}
}

void MainWindow::openProcessingSystemTab(QString system_name)
{
	NGSD db;
	int sys_id = db.processingSystemId(system_name, false);
	if (sys_id==-1)
	{
		GUIHelper::showMessage("NGSD error", "Processing system name '" + system_name + "' not found in NGSD!");
		return;
	}

	ProcessingSystemWidget* widget = new ProcessingSystemWidget(this, sys_id);
	int index = openTab(QIcon(":/Icons/NGSD_processing_system.png"), db.getProcessingSystemData(sys_id).name, widget);
	if (Settings::boolean("debug_mode_enabled"))
	{
		ui_.tabs->setTabToolTip(index, "NGSD ID: " + QString::number(sys_id));
	}
}

void MainWindow::openProjectTab(QString name)
{
	ProjectWidget* widget = new ProjectWidget(this, name);
	int index = openTab(QIcon(":/Icons/NGSD_project.png"), name, widget);
	if (Settings::boolean("debug_mode_enabled"))
	{
		ui_.tabs->setTabToolTip(index, "NGSD ID: " + NGSD().getValue("SELECT id FROM project WHERE name=:0", true, name).toString());
	}
}

int MainWindow::openTab(QIcon icon, QString name, QWidget* widget)
{
	QScrollArea* scroll_area = new QScrollArea(this);
	scroll_area->setFrameStyle(QFrame::NoFrame);
	scroll_area->setWidgetResizable(true);
	scroll_area->setWidget(widget);
	//fix color problems
	QPalette pal;
	pal.setColor(QPalette::Window,QColor(0,0,0,0));
	scroll_area->setPalette(pal);
	scroll_area->setBackgroundRole(QPalette::Window);
	scroll_area->widget()->setPalette(pal);
	scroll_area->widget()->setBackgroundRole(QPalette::Window);
	//show tab
	int index = ui_.tabs->addTab(scroll_area, icon, name);
	ui_.tabs->setCurrentIndex(index);

	return index;
}

void MainWindow::closeTab(int index)
{
	if (index==0)
	{
		int res = QMessageBox::question(this, "Close file?", "Do you want to close the current sample?", QMessageBox::Yes, QMessageBox::No);
		if (res==QMessageBox::Yes)
		{
			loadFile();
		}
	}
	else
	{
		QWidget* widget = ui_.tabs->widget(index);
		ui_.tabs->removeTab(index);
		widget->deleteLater();
	}
}

void MainWindow::on_actionChangeLog_triggered()
{
	QDesktopServices::openUrl(QUrl("https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/changelog.md"));
}

void MainWindow::loadFile(QString filename, bool show_only_error_issues)
{
	//store variant list in case it changed
	if (!variants_changed_.isEmpty())
	{
		int result = QMessageBox::question(this, "Store GSvar file?", "The GSvar file was changed by you.\nDo you want to store the changes to file?", QMessageBox::Yes, QMessageBox::No);
		if (result==QMessageBox::Yes)
		{
			storeCurrentVariantList();
		}
	}

	QTime timer;
	timer.start();

	//reset GUI and data structures
	setWindowTitle(appName());
	filename_ = "";
	variants_.clear();
	GlobalServiceProvider::clearFileLocationProvider();
	variants_changed_.clear();
	cnvs_.clear();
	svs_.clear();
	IgvSessionManager::clearAll();
	ui_.vars->clearContents();
	report_settings_ = ReportSettings();
	connect(report_settings_.report_config.data(), SIGNAL(variantsChanged()), this, SLOT(storeReportConfig()));
	germline_report_ps_ = "";
	somatic_report_settings_ = SomaticReportSettings();

	ui_.tabs->setCurrentIndex(0);

	ui_.filters->reset(true);

	Log::perf("Clearing variant table took ", timer);

	if (filename=="") return;

	//load data
	QApplication::setOverrideCursor(Qt::BusyCursor);
	try
	{
		//load variants
		timer.restart();
		variants_.load(filename);

		Log::perf("Loading small variant list took ", timer);
		QString mode_title = "";
		if (Helper::isHttpUrl(filename))
		{
			GlobalServiceProvider::setFileLocationProvider(QSharedPointer<FileLocationProviderRemote>(new FileLocationProviderRemote(filename)));
		}
		else
		{
			GlobalServiceProvider::setFileLocationProvider(QSharedPointer<FileLocationProviderLocal>(new FileLocationProviderLocal(filename, variants_.getSampleHeader(), variants_.type())));
			mode_title = " (local mode)";
		}

		//load CNVs
		timer.restart();
		FileLocation cnv_loc = GlobalServiceProvider::fileLocationProvider().getAnalysisCnvFile();
		if (cnv_loc.exists)
		{
			try
			{
				cnvs_.load(cnv_loc.filename);
				int cnv_count_initial = cnvs_.count();
				double min_ll = 0.0;
				while (cnvs_.count()>50000)
				{
					min_ll += 1.0;
					FilterResult result(cnvs_.count());
					FilterCnvLoglikelihood filter;
					filter.setDouble("min_ll", min_ll);
					filter.apply(cnvs_, result);
					result.removeFlagged(cnvs_);
				}
				if (min_ll>0)
				{
					QMessageBox::information(this, "CNV pre-filtering applied", "The CNV calls file contains too many CNVs: " + QString::number(cnv_count_initial) +".\nOnly CNVs with log-likelyhood >= " + QString::number(min_ll) +" are displayed: " + QString::number(cnvs_.count()) +".");
				}
			}
			catch(Exception& e)
			{
				QMessageBox::warning(this, "Error loading CNVs", e.message());
				cnvs_.clear();
			}
		}
		Log::perf("Loading CNV list took ", timer);

		//load SVs
		timer.restart();
		FileLocation sv_loc = GlobalServiceProvider::fileLocationProvider().getAnalysisSvFile();
		if (sv_loc.exists)
		{
			try
			{
				svs_.load(sv_loc.filename);
			}
			catch(Exception& e)
			{
				QMessageBox::warning(this, "Error loading SVs", e.message());
				svs_.clear();
			}
		}
		Log::perf("Loading SV list took ", timer);

		//determine valid filter entries from filter column (and add new filters low_mappability/mosaic to make outdated GSvar files work as well)
		QStringList valid_filter_entries = variants_.filters().keys();
		if (!valid_filter_entries.contains("low_mappability")) valid_filter_entries << "low_mappability";
		if (!valid_filter_entries.contains("mosaic")) valid_filter_entries << "mosaic";
		ui_.filters->setValidFilterEntries(valid_filter_entries);

		//update data structures
		Settings::setPath("path_variantlists", filename);
		filename_ = filename;

		//update GUI
		setWindowTitle(appName() + " - " + variants_.analysisName() + mode_title);
		ui_.statusBar->showMessage("Loaded variant list with " + QString::number(variants_.count()) + " variants.");

		refreshVariantTable(false);
		ui_.vars->adaptColumnWidths();

		QApplication::restoreOverrideCursor();
	}
	catch(Exception& e)
	{
		QApplication::restoreOverrideCursor();
		QMessageBox::warning(this, "Error", "Loading the file '" + filename + "' or displaying the contained variants failed!\nError message:\n" + e.message());
		loadFile();
		return;
	}

	//check analysis for issues (outdated, missing columns, wrong genome build, bad quality, ...)
	QList<QPair<Log::LogLevel, QString>> issues;
	try
	{
		ui_.variant_details->setLabelTooltips(variants_);
	}
	catch(Exception& e)
	{
		issues << qMakePair(Log::LOG_INFO, e.message());
	}
	checkVariantList(issues);

	//check variant list in NGSD
	checkProcessedSamplesInNGSD(issues);

	//show issues
	if (showAnalysisIssues(issues, show_only_error_issues)==QDialog::Rejected)
	{
		loadFile();
		return;
	}

	//load report config
	if (germlineReportSupported())
	{
		loadReportConfig();
	}
	else if(LoginManager::active() && somaticReportSupported())
	{
		loadSomaticReportConfig();
	}

	//check mendelian error rate for trios
	AnalysisType type = variants_.type();
	if (type==GERMLINE_TRIO)
	{
		checkMendelianErrorRate();
	}

	//notify for variant validation
	checkPendingVariantValidations();

	//activate Circos plot menu item if plot is available
	if (type==GERMLINE_SINGLESAMPLE && !GlobalServiceProvider::fileLocationProvider().getCircosPlotFiles(false).isEmpty())
	{
		ui_.actionCircos->setEnabled(true);
	}
	else
	{
		ui_.actionCircos->setEnabled(false);
	}

	//activate Repeat Expansion menu item if RE calls are available
	if (type==GERMLINE_SINGLESAMPLE && !GlobalServiceProvider::fileLocationProvider().getRepeatExpansionFiles(false).isEmpty())
	{
		ui_.actionRE->setEnabled(true);
	}
	else
	{
		ui_.actionRE->setEnabled(false);
	}

	//activate PRS menu item if PRS are available
	if (type==GERMLINE_SINGLESAMPLE && !GlobalServiceProvider::fileLocationProvider().getPrsFiles(false).isEmpty())
	{
		ui_.actionPRS->setEnabled(true);
	}
	else
	{
		ui_.actionPRS->setEnabled(false);
	}

	//activate virus table
	ui_.actionVirusDetection->setEnabled(false);
	if (type==SOMATIC_PAIR && NGSD::isAvailable())
	{
		QString ps_tumor = variants_.mainSampleName();
		NGSD db;
		QString ps_tumor_id = db.processedSampleId(ps_tumor, false);
		if (GlobalServiceProvider::database().processedSamplePath(ps_tumor_id, PathType::VIRAL).exists)
		{
			ui_.actionVirusDetection->setEnabled(true);
		}
	}

	//activate cfDNA menu entries and get all available cfDNA samples
	cf_dna_available = false;
	ui_.actionDesignCfDNAPanel->setVisible(false);
	ui_.actionCfDNADiseaseCourse->setVisible(false);
	ui_.actionDesignCfDNAPanel->setEnabled(false);
	ui_.actionCfDNADiseaseCourse->setEnabled(false);
	cfdna_menu_btn_->setVisible(false);
	cfdna_menu_btn_->setEnabled(false);
	if (somaticReportSupported() || tumoronlyReportSupported())
	{
		ui_.actionDesignCfDNAPanel->setVisible(true);
		ui_.actionCfDNADiseaseCourse->setVisible(true);
		cfdna_menu_btn_->setVisible(true);

		if (LoginManager::active())
		{
			ui_.actionDesignCfDNAPanel->setEnabled(true);
			cfdna_menu_btn_->setEnabled(true);
			NGSD db;

			// get all same samples
			int sample_id = db.sampleId(variants_.mainSampleName()).toInt();
			QSet<int> same_sample_ids = db.relatedSamples(sample_id, "same sample");
			same_sample_ids << sample_id; // add current sample id

			// get all related cfDNA
			QSet<int> cf_dna_sample_ids;
			foreach (int cur_sample_id, same_sample_ids)
			{
				cf_dna_sample_ids.unite(db.relatedSamples(cur_sample_id, "tumor-cfDNA"));
			}

			if (cf_dna_sample_ids.size() > 0)
			{
				ui_.actionCfDNADiseaseCourse->setEnabled(true);
				cf_dna_available = true;
			}
		}
	}

	//activate RNA menu
	ui_.actionExpressionData->setEnabled(false);
	ui_.actionExonExpressionData->setEnabled(false);
	ui_.actionShowSplicing->setEnabled(false);
	ui_.actionShowRnaFusions->setEnabled(false);
	if (LoginManager::active())
	{
		NGSD db;

		QString sample_id = db.sampleId(filename_, false);
		if (sample_id!="")
		{
			foreach (int rna_sample_id, db.relatedSamples(sample_id.toInt(), "same sample", "RNA"))
			{
				// check for required files
				foreach (const QString& rna_ps_id, db.getValues("SELECT id FROM processed_sample WHERE sample_id=:0", QString::number(rna_sample_id)))
				{
					if (!db.userCanAccess(LoginManager::userId(), rna_ps_id.toInt())) continue;

					// search for gene count file
					FileLocation rna_count_file = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::EXPRESSION);
					if (rna_count_file.exists) ui_.actionExpressionData->setEnabled(true);

					// search for gene count file
					rna_count_file = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::EXPRESSION_EXON);
					if (rna_count_file.exists) ui_.actionExonExpressionData->setEnabled(true);

					// search for splicing predictions tsv
					FileLocation splicing_file = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::SPLICING_ANN);
					if (splicing_file.exists) ui_.actionShowSplicing->setEnabled(true);


					// search for arriba fusion file
					FileLocation arriba_fusion_file = GlobalServiceProvider::database().processedSamplePath(rna_ps_id, PathType::FUSIONS);
					if (arriba_fusion_file.exists) ui_.actionShowRnaFusions->setEnabled(true);
				}
			}
		}
	}
}

void MainWindow::checkVariantList(QList<QPair<Log::LogLevel, QString>>& issues)
{
	//check genome builds match
	if (variants_.build()!=GSvarHelper::build())
	{
		issues << qMakePair(Log::LOG_ERROR, "Genome build of GSvar file (" + buildToString(variants_.build(), true) + ") not matching genome build of the GSvar application (" + buildToString(GSvarHelper::build(), true) + ")! Re-do the analysis for " + buildToString(GSvarHelper::build(), true) +"!");
	}

	//check creation date
	QDate create_date = variants_.getCreationDate();
	if (create_date.isValid())
	{
		if (create_date < QDate::currentDate().addDays(-42))
		{
			issues << qMakePair(Log::LOG_INFO, "Variant annotations are older than six weeks (" + create_date.toString("yyyy-MM-dd") + ").");
		}
		QDate gsvar_file_outdated_before = QDate::fromString(Settings::string("gsvar_file_outdated_before", true), "yyyy-MM-dd");
		if (gsvar_file_outdated_before.isValid() && create_date<gsvar_file_outdated_before)
		{
			issues << qMakePair(Log::LOG_WARNING, "Variant annotations are outdated! They are older than " + gsvar_file_outdated_before.toString("yyyy-MM-dd") + ". Please re-annotate variants!");
		}
	}

	//check sample header
	try
	{
		variants_.getSampleHeader();
	}
	catch(Exception e)
	{
		issues << qMakePair(Log::LOG_WARNING,  e.message());
	}

	//create list of required columns
	QStringList cols;
	cols << "filter";
	AnalysisType type = variants_.type();
	if (type==GERMLINE_SINGLESAMPLE || type==GERMLINE_TRIO || type==GERMLINE_MULTISAMPLE)
	{
		cols << "classification";
		cols << "NGSD_hom";
		cols << "NGSD_het";
		cols << "comment";
		cols << "gene_info";
	}
	if (type==SOMATIC_SINGLESAMPLE || type==SOMATIC_PAIR || type==CFDNA)
	{
		cols << "somatic_classification";
		cols << "somatic_classification_comment";
		cols << "NGSD_som_vicc_interpretation";
		cols << "NGSD_som_vicc_comment";
	}

	//check columns
	foreach(const QString& col, cols)
	{
		if (variants_.annotationIndexByName(col, true, false)==-1)
		{
			issues << qMakePair(Log::LOG_WARNING, "Column '" + col + "' missing. Please re-annotate variants!");
		}
	}

	//check data was loaded completely
	if (germlineReportSupported())
	{
		NGSD db;
		int sys_id = db.processingSystemIdFromProcessedSample(germlineReportSample());
		ProcessingSystemData sys_data = db.getProcessingSystemData(sys_id);
		if (sys_data.type=="WES" || sys_data.type=="WGS" || sys_data.type=="lrGS")
		{
			QSet<Chromosome> chromosomes;
			for(int i=0; i<variants_.count(); ++i)
			{
				chromosomes << variants_[i].chr();
			}
			if (chromosomes.size()<23)
			{
				issues << qMakePair(Log::LOG_WARNING, "Variants detected on " + QString::number(chromosomes.size()) + " chromosomes only! Expected variants on at least 23 chromosomes for WES/WGS data! Please re-do variant calling of small variants!");
			}
		}
	}

	//check sv annotation
	if (type==GERMLINE_SINGLESAMPLE || type==GERMLINE_TRIO || type==GERMLINE_MULTISAMPLE)
	{
		if (svs_.isValid())
		{
			//check for NGSD count annotation
			if(!svs_.annotationHeaders().contains("NGSD_HOM") || !svs_.annotationHeaders().contains("NGSD_HET") || !svs_.annotationHeaders().contains("NGSD_AF"))
			{
				issues << qMakePair(Log::LOG_WARNING, QString("Annotation of structural variants is outdated! Please re-annotate structural variants!"));
			}
		}
	}

	//check GenLab Labornummer is present (for diagnostic samples only)
	if (GenLabDB::isAvailable() && NGSD::isAvailable())
	{
		NGSD db;
		foreach(const SampleInfo& info, variants_.getSampleHeader())
		{
			QString ps_name = info.id;
			QString ps_id = db.processedSampleId(ps_name, false);
			if (ps_id=="") continue; //not in NGSD

			QString project_type = db.getValue("SELECT p.type FROM processed_sample ps, project p WHERE p.id=ps.project_id AND ps.id=" + ps_id).toString();
			if (project_type=="diagnostic")
			{
				GenLabDB genlab;
				QString genlab_patient_id = genlab.patientIdentifier(ps_name);
				if (genlab_patient_id.isEmpty())
				{
					issues << qMakePair(Log::LOG_WARNING, QString("GenLab Labornummer probably not set correctly for dianostic sample '" + ps_name + "'. Please correct the Labornummer in GenLab!"));
				}
			}
		}
	}
}

void MainWindow::checkProcessedSamplesInNGSD(QList<QPair<Log::LogLevel, QString>>& issues)
{
	if (!LoginManager::active()) return;

	NGSD db;

	foreach(const SampleInfo& info, variants_.getSampleHeader())
	{
		QString ps = info.id;
		QString ps_id = db.processedSampleId(ps, false);
		if (ps_id=="") continue;

		//check quality
		QString quality = db.getValue("SELECT quality FROM processed_sample WHERE id=" + ps_id).toString();
		if (quality=="bad")
		{
			issues << qMakePair(Log::LOG_WARNING, "Quality of processed sample '" + ps + "' is 'bad'!");
		}
		else if (quality=="n/a")
		{
			issues << qMakePair(Log::LOG_WARNING, "Quality of processed sample '" + ps + "' is not set!");
		}

		//check sequencing run is marked as analyzed
		QString run_status = db.getValue("SELECT r.status FROM sequencing_run r, processed_sample ps WHERE r.id=ps.sequencing_run_id AND ps.id=" + ps_id).toString();
		if (run_status!="analysis_finished")
		{
			issues << qMakePair(Log::LOG_WARNING, "Sequencing run of the sample '" + ps + "' does not have status 'analysis_finished'!");
		}

		//check KASP result
		try
		{
			double error_prob = db.kaspData(ps_id).random_error_prob;
			if (error_prob>0.03)
			{
				issues << qMakePair(Log::LOG_WARNING, "KASP swap probability of processed sample '" + ps + "' is larger than 3%!");
			}
		}
		catch (DatabaseException /*e*/)
		{
			//nothing to do (KASP not done or invalid)
		}

		//check variants are imported
		AnalysisType type = variants_.type();
		if (type==GERMLINE_SINGLESAMPLE || type==GERMLINE_TRIO || type==GERMLINE_MULTISAMPLE)
		{

			QString sys_type = db.getProcessingSystemData(db.processingSystemIdFromProcessedSample(ps)).type;
			if (sys_type=="WGS" || sys_type=="WES" || sys_type=="lrGS")
			{
				ImportStatusGermline import_status = db.importStatus(ps_id);
				if (import_status.small_variants==0)
				{
					issues << qMakePair(Log::LOG_WARNING, "No germline variants imported into NGSD for processed sample '" + ps + "'!");
				}
			}
		}
	}
}

int MainWindow::showAnalysisIssues(QList<QPair<Log::LogLevel, QString>>& issues, bool show_only_error_issues)
{
	if (issues.empty()) return QDialog::Accepted;

	//generate text
	bool has_error = false;
	QStringList lines;
	foreach(auto issue, issues)
	{
		if (issue.first==Log::LOG_ERROR)
		{
			has_error = true;
			lines << "<font color=red>Error:</font>";
		}
		else if (issue.first==Log::LOG_WARNING)
		{
			lines << "<font color=orange>Warning:</font>";
		}
		else
		{
			lines << "Notice:";
		}
		lines << issue.second;
		lines << "";
	}

	if (show_only_error_issues && !has_error) return QDialog::Accepted;

	//show dialog
	QLabel* label = new QLabel(lines.join("<br>"));
	label->setMargin(6);
	auto dlg = GUIHelper::createDialog(label, "GSvar analysis issues", "The following issues were encountered when loading the analysis:", true);
	return dlg->exec();
}

void MainWindow::on_actionAbout_triggered()
{
	QString about_text = appName()+ " " + QCoreApplication::applicationVersion();

	about_text += "\n\nA free viewing and filtering tool for genomic variants.";

	//general infos
	about_text += "\n";
	about_text += "\nGenome build: " + buildToString(GSvarHelper::build());
	about_text += "\nArchitecture: " + QSysInfo::buildCpuArchitecture();

	//client-server infos
	about_text += "\n";
	if (ClientHelper::isClientServerMode())
	{
		about_text += "\nMode: client-server";
        int status_code = -1;
        ServerInfo server_info = ClientHelper::getServerInfo(status_code);
        if (status_code!=200)
        {
            about_text += "\nServer returned " + QString::number(status_code);
        }
        else
        {
            about_text += "\nServer version: " + server_info.version;
            about_text += "\nAPI version: " + server_info.api_version;
            about_text += "\nServer start time: " + server_info.server_start_time.toString("yyyy-MM-dd hh:mm:ss");
        }
    }
	else
	{
		about_text += "\nMode: stand-alone (no server)";
	}

	QMessageBox::about(this, "About " + appName(), about_text);
}

void MainWindow::loadReportConfig()
{
	//check if applicable
	if (!germlineReportSupported()) return;

	//check if report config exists
	NGSD db;
	QString processed_sample_id = db.processedSampleId(germlineReportSample(), false);
	int rc_id = db.reportConfigId(processed_sample_id);
	if (rc_id==-1) return;

	//load
	report_settings_.report_config = db.reportConfig(rc_id, variants_, cnvs_, svs_);
	connect(report_settings_.report_config.data(), SIGNAL(variantsChanged()), this, SLOT(storeReportConfig()));


	//updateGUI
	refreshVariantTable();
}

void MainWindow::loadSomaticReportConfig()
{
	if(filename_ == "") return;

	NGSD db;

	//Determine processed sample ids
	QString ps_tumor = variants_.mainSampleName();
	QString ps_tumor_id = db.processedSampleId(ps_tumor, false);
	if(ps_tumor_id == "") return;
	QString ps_normal = normalSampleName();
	if(ps_normal.isEmpty()) return;
	QString ps_normal_id = db.processedSampleId(ps_normal, false);
	if(ps_normal_id == "") return;

	somatic_report_settings_.tumor_ps = ps_tumor;
	somatic_report_settings_.normal_ps = ps_normal;
	somatic_report_settings_.msi_file = GlobalServiceProvider::fileLocationProvider().getSomaticMsiFile().filename;
	somatic_report_settings_.viral_file = GlobalServiceProvider::database().processedSamplePath(ps_tumor_id, PathType::VIRAL).filename;

	somatic_report_settings_.sbs_signature = GlobalServiceProvider::fileLocationProvider().getSignatureSbsFile().filename;
	somatic_report_settings_.id_signature = GlobalServiceProvider::fileLocationProvider().getSignatureIdFile().filename;
	somatic_report_settings_.dbs_signature = GlobalServiceProvider::fileLocationProvider().getSignatureDbsFile().filename;
	somatic_report_settings_.cnv_signature = GlobalServiceProvider::fileLocationProvider().getSignatureCnvFile().filename;

	try //load normal sample
	{
		somatic_control_tissue_variants_.load(GlobalServiceProvider::database().processedSamplePath(db.processedSampleId(ps_normal), PathType::GSVAR).filename);
	}
	catch(Exception e)
	{
		QMessageBox::warning(this, "Could not load germline GSvar file", "Could not load germline GSvar file. No germline variants will be parsed for somatic report generation. Message: " + e.message());
	}



	//Continue loading report (only if existing in NGSD)
	if(db.somaticReportConfigId(ps_tumor_id, ps_normal_id) == -1) return;


	QStringList messages;
	somatic_report_settings_.report_config = db.somaticReportConfig(ps_tumor_id, ps_normal_id, variants_, cnvs_, somatic_control_tissue_variants_, messages);


	if(!messages.isEmpty())
	{
		QMessageBox::warning(this, "Somatic report configuration", "The following problems were encountered while loading the som. report configuration:\n" + messages.join("\n"));
	}

	//Preselect target region bed file in NGSD
	if(somatic_report_settings_.report_config.targetRegionName()!="")
	{
		ui_.filters->setTargetRegionByDisplayName(somatic_report_settings_.report_config.targetRegionName());
	}

	//Preselect filter from NGSD som. rep. conf.
	if(somatic_report_settings_.report_config.filter() != "")
	{
		ui_.filters->setFilter( somatic_report_settings_.report_config.filter() );
	}

	somatic_report_settings_.target_region_filter = ui_.filters->targetRegion();

	refreshVariantTable();
}

void MainWindow::storeSomaticReportConfig()
{
	if(filename_ == "") return;
	if(!LoginManager::active()) return;
	if(variants_.type() != SOMATIC_PAIR) return;

	NGSD db;
	QString ps_tumor_id = db.processedSampleId(variants_.mainSampleName(), false);
	QString ps_normal_id = db.processedSampleId(normalSampleName(), false);

	if(ps_tumor_id=="" || ps_normal_id == "")
	{
		QMessageBox::warning(this, "Storing somatic report configuration", "Samples were not found in the NGSD!");
		return;
	}

	int conf_id = db.somaticReportConfigId(ps_tumor_id, ps_normal_id);

	if (conf_id!=-1)
	{
		SomaticReportConfigurationData conf_creation = db.somaticReportConfigData(conf_id);
		if (conf_creation.last_edit_by!="" && conf_creation.last_edit_by!=LoginManager::userName())
		if (QMessageBox::question(this, "Storing report configuration", conf_creation.history() + "\n\nDo you want to update/override it?")==QMessageBox::No)
		{
			return;
		}
	}

	//store
	try
	{
		db.setSomaticReportConfig(ps_tumor_id, ps_normal_id, somatic_report_settings_.report_config, variants_, cnvs_, somatic_control_tissue_variants_, Helper::userName());
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Storing somatic report configuration", "Error: Could not store the somatic report configuration.\nPlease resolve this error or report it to the administrator:\n\n" + e.message());
	}
}

void MainWindow::storeReportConfig()
{
	//check if applicable
	if (!germlineReportSupported()) return;

	//check sample
	NGSD db;
	QString processed_sample_id = db.processedSampleId(germlineReportSample(), false);
	if (processed_sample_id=="")
	{
		QMessageBox::warning(this, "Storing report configuration", "Sample was not found in the NGSD!");
		return;
	}

	//check if config exists and not edited by other user
	int conf_id = db.reportConfigId(processed_sample_id);
	if (conf_id!=-1)
	{
		QSharedPointer<ReportConfiguration> report_config = db.reportConfig(conf_id, variants_, cnvs_, svs_);
		if (report_config->lastUpdatedBy()!="" && report_config->lastUpdatedBy()!=LoginManager::userName())
		{
			if (QMessageBox::question(this, "Storing report configuration", report_config->history() + "\n\nDo you want to override it?")==QMessageBox::No)
			{
				return;
			}
		}
	}

	//store
	try
	{
		report_settings_.report_config.data()->blockSignals(true); //block signals - otherwise the variantsChanged signal is emitted and storeReportConfig is called again, which leads to hanging of the application because of database locks
		db.setReportConfig(processed_sample_id, report_settings_.report_config, variants_, cnvs_, svs_);
		report_settings_.report_config.data()->blockSignals(false);
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Storing report configuration", e.message());
	}
}

void MainWindow::generateEvaluationSheet()
{
	//check if applicable
	if (!germlineReportSupported()) return;

	QString base_name = germlineReportSample();

	//make sure free-text phenotype infos are available
	NGSD db;
	QString sample_id = db.sampleId(base_name);
	QList<SampleDiseaseInfo> disease_infos = db.getSampleDiseaseInfo(sample_id, "clinical phenotype (free text)");
	if (disease_infos.isEmpty() && QMessageBox::question(this, "Evaluation sheet", "No clinical phenotype (free text) is set for the sample!\nIt will be shown on the evaluation sheet!\n\nDo you want to set it?")==QMessageBox::Yes)
	{
		SampleDiseaseInfoWidget* widget = new SampleDiseaseInfoWidget(base_name, this);
		widget->setDiseaseInfo(db.getSampleDiseaseInfo(sample_id));
		auto dlg = GUIHelper::createDialog(widget, "Sample disease details", "", true);
		if (dlg->exec() != QDialog::Accepted) return;

		db.setSampleDiseaseInfo(sample_id, widget->diseaseInfo());
	}

	//try to get VariantListInfo from the NGSD
	QString ps_id = db.processedSampleId(base_name);
	EvaluationSheetData evaluation_sheet_data = db.evaluationSheetData(ps_id, false);
	evaluation_sheet_data.build = GSvarHelper::build();
	if (evaluation_sheet_data.ps_id == "") //No db entry found > init
	{
		evaluation_sheet_data.ps_id = db.processedSampleId(base_name);
		evaluation_sheet_data.dna_rna = db.getSampleData(sample_id).name_external;
		// make sure reviewer 1 contains name not user id
		evaluation_sheet_data.reviewer1 = db.userName(db.userId(report_settings_.report_config->createdBy()));
		evaluation_sheet_data.review_date1 = report_settings_.report_config->createdAt().date();
		evaluation_sheet_data.reviewer2 = LoginManager::userName();
		evaluation_sheet_data.review_date2 = QDate::currentDate();
	}

	//Show VariantSheetEditDialog
	EvaluationSheetEditDialog* edit_dialog = new EvaluationSheetEditDialog(this);
	edit_dialog->importEvaluationSheetData(evaluation_sheet_data);
	if (edit_dialog->exec() != QDialog::Accepted) return;

	//Store updated info in the NGSD
	db.storeEvaluationSheetData(evaluation_sheet_data, true);

	//get filename
	QString folder = Settings::path("gsvar_variantsheet_folder");
	QString filename = QFileDialog::getSaveFileName(this, "Store evaluation sheet",  folder + "/" + base_name + "_variant_sheet_" + QDate::currentDate().toString("yyyyMMdd") + ".html", "HTML files (*.html);;All files(*.*)");
	if (filename.isEmpty()) return;

	//write sheet
	PrsTable prs_table; //not needed
	GermlineReportGeneratorData generator_data(GSvarHelper::build(), base_name, variants_, cnvs_, svs_, prs_table, report_settings_, ui_.filters->filters(), GSvarHelper::preferredTranscripts(), GlobalServiceProvider::statistics());
	GermlineReportGenerator generator(generator_data);
	generator.writeEvaluationSheet(filename, evaluation_sheet_data);

	if (QMessageBox::question(this, "Evaluation sheet", "Evaluation sheet generated successfully!\nDo you want to open it in your browser?")==QMessageBox::Yes)
	{
		QDesktopServices::openUrl(QUrl::fromLocalFile(filename));
	}
}

void MainWindow::transferSomaticData()
{
	try
	{
		if(variants_.type()!=AnalysisType::SOMATIC_PAIR)
		{
			INFO(Exception, "Error: only possible for tumor-normal pair!");
		}

		SomaticDataTransferWidget data_transfer(somatic_report_settings_.tumor_ps, somatic_report_settings_.normal_ps, this);
		data_transfer.exec();
	}
	catch(Exception e)
	{
		GUIHelper::showException(this, e, "Transfer somatic data to MTB");
	}
}

void MainWindow::showReportConfigInfo()
{
	//check if applicable
	if (!germlineReportSupported() && !somaticReportSupported()) return;

	QString ps = germlineReportSupported() ? germlineReportSample() : variants_.mainSampleName();
	QString title = "Report configuration information of " + ps;

	//check sample exists
	NGSD db;
	QString processed_sample_id = db.processedSampleId(ps, false);
	if (processed_sample_id=="")
	{
		QMessageBox::warning(this, title, "Sample was not found in the NGSD!");
		return;
	}

	//check config exists
	if(germlineReportSupported())
	{
		int conf_id = db.reportConfigId(processed_sample_id);
		if (conf_id==-1)
		{
			QMessageBox::warning(this, title , "No germline report configuration found in the NGSD!");
			return;
		}


		QSharedPointer<ReportConfiguration> report_config = db.reportConfig(conf_id, variants_, cnvs_, svs_);

		QMessageBox::information(this, title, report_config->history() + "\n\n" + report_config->variantSummary());
	}
	else if(somaticReportSupported())
	{
		QString ps_normal_id = db.processedSampleId(normalSampleName(), false);
		int conf_id = db.somaticReportConfigId(processed_sample_id, ps_normal_id);
		if(conf_id==-1)
		{
			QMessageBox::warning(this, title, "No somatic report configuration found in the NGSD!");
			return;
		}
		QMessageBox::information(this, title, db.somaticReportConfigData(conf_id).history());
	}
}

void MainWindow::editOtherCausalVariant()
{
	QString title = "Add/edit other causal variant";
	try
	{
		//check if applicable
		if (!germlineReportSupported()) INFO(ArgumentException, "This feature is only available for germline!");

		QString ps = germlineReportSample();

		//check sample exists in NGSD
		NGSD db;
		QString processed_sample_id = db.processedSampleId(ps, false);
		if (processed_sample_id=="") INFO(ArgumentException, "Sample was not found in the NGSD!");

		// get report config
		OtherCausalVariant causal_variant = report_settings_.report_config->otherCausalVariant();
		QStringList variant_types = db.getEnum("report_configuration_other_causal_variant", "type");
		QStringList inheritance_modes = db.getEnum("report_configuration_other_causal_variant", "inheritance");

		//open edit dialog
		CausalVariantEditDialog dlg(causal_variant, variant_types, inheritance_modes, this);
		dlg.setWindowTitle(title + " of " + ps);

		if (dlg.exec()!=QDialog::Accepted) return;

		//store updated causal variant in NGSD
		if (dlg.causalVariant().isValid())
		{
			report_settings_.report_config->setOtherCausalVariant(dlg.causalVariant());
			report_settings_.report_config->variantsChanged();
		}
	}
	catch(Exception e)
	{
		GUIHelper::showException(this, e, title);
	}
}

void MainWindow::deleteOtherCausalVariant()
{
	QString title = "Delete other causal variant";
	try
	{
		//check if applicable
		if (!germlineReportSupported()) INFO(ArgumentException, "This feature is only available for germline!");

		QString ps = germlineReportSample();

		//check sample exists
		NGSD db;
		QString processed_sample_id = db.processedSampleId(ps, false);
		if (processed_sample_id=="") INFO(ArgumentException, "Sample was not found in the NGSD!");

		OtherCausalVariant causal_variant = report_settings_.report_config->otherCausalVariant();
		if(!causal_variant.isValid()) return;

		//show dialog to confirm by user
		QString message_text = "Are you sure you want to delete the following causal variant?\n" + causal_variant.type + " at " + causal_variant.coordinates + " (gene: " + causal_variant.gene + ", comment: " + causal_variant.comment.replace("\n", " ") + ")";
		QMessageBox::StandardButton response = QMessageBox::question(this, title + " of " + ps, message_text, QMessageBox::Yes|QMessageBox::No);
		if(response != QMessageBox::Yes) return;

		//replace other causal variant with empty struct
		report_settings_.report_config->setOtherCausalVariant(OtherCausalVariant());
	}
	catch(Exception e)
	{
		GUIHelper::showException(this, e, title);
	}
}

void MainWindow::finalizeReportConfig()
{
	QString title = "Finalize report configuration";
	try
	{
		//check if applicable
		if(!germlineReportSupported()) INFO(ArgumentException, "This feature is only available for germline!");

		//check sample exists
		NGSD db;
		QString processed_sample_id = db.processedSampleId(germlineReportSample(), false);
		if (processed_sample_id=="") INFO(ArgumentException, "Sample was not found in the NGSD!");

		//check config exists
		int conf_id = db.reportConfigId(processed_sample_id);
		if (conf_id==-1) INFO(ArgumentException, "No report configuration for this sample found in the NGSD!");

		//make sure the user knows what he does
		int button = QMessageBox::question(this, title, "Do you really want to finalize the report configuration?\nIt cannot be modified or deleted when finalized!");
		if (button!=QMessageBox::Yes) return;

		//finalize
		db.finalizeReportConfig(conf_id, LoginManager::userId());

		//update report settings data structure
		report_settings_.report_config = db.reportConfig(conf_id, variants_, cnvs_, svs_);
		connect(report_settings_.report_config.data(), SIGNAL(variantsChanged()), this, SLOT(storeReportConfig()));
	}
	catch(Exception e)
	{
		GUIHelper::showException(this, e, title);
	}
}

void MainWindow::generateReport()
{
	if (filename_=="") return;

	//check if this is a germline or somatic
	if (somaticReportSupported())
	{
		generateReportSomaticRTF();
	}
	else if (tumoronlyReportSupported())
	{
		generateReportTumorOnly();
	}
	else if (germlineReportSupported())
	{
		generateReportGermline();
	}
	else
	{
		QMessageBox::information(this, "Report error", "Report not supported for this type of analysis!");
	}
}


void MainWindow::generateReportTumorOnly()
{
	try
	{
		TumorOnlyReportWorker::checkAnnotation(variants_);
	}
	catch(FileParseException e)
	{
		QMessageBox::warning(this, "Invalid tumor only file" + filename_, "Could not find all neccessary annotations in tumor-only variant list. Aborting creation of report. " + e.message());
		return;
	}
	QString ps = variants_.mainSampleName();

	NGSD db;

	//get report settings
	TumorOnlyReportWorkerConfig config;
	config.threads = Settings::integer("threads");
	int sys_id = db.processingSystemIdFromProcessedSample(ps);

	config.sys = db.getProcessingSystemData(sys_id);
	config.ps_data = db.getProcessedSampleData(db.processedSampleId(ps));
	config.roi = ui_.filters->targetRegion();

	config.low_coverage_file = GlobalServiceProvider::fileLocationProvider().getSomaticLowCoverageFile().filename;
	config.bam_file = GlobalServiceProvider::fileLocationProvider().getBamFiles(true).at(0).filename;
	config.filter_result = filter_result_;
	config.preferred_transcripts = GSvarHelper::preferredTranscripts();
	config.build = GSvarHelper::build();

	TumorOnlyReportDialog dlg(variants_, config, this);
	if(!dlg.exec()) return;

	//get RTF file name
	QString destination_path = last_report_path_ + "/" + ps + "_DNA_tumor_only_" + QDate::currentDate().toString("yyyyMMdd") + ".rtf";
	QString file_rep = QFileDialog::getSaveFileName(this, "Store report file", destination_path, "RTF files (*.rtf);;All files(*.*)");
	if (file_rep=="") return;

	//Generate RTF
	QApplication::setOverrideCursor(Qt::BusyCursor);
	try
	{
		TumorOnlyReportWorker worker(variants_, config);

		QByteArray temp_filename = Helper::tempFileName(".rtf").toUtf8();
		worker.writeRtf(temp_filename);
		ReportWorker::moveReport(temp_filename, file_rep);

		if(!ui_.filters->targetRegion().isValid()) //if no ROI filter was set, use panel target information instead
		{
			TargetRegionInfo roi_info;
			roi_info.name = config.sys.name;
			roi_info.regions = GlobalServiceProvider::database().processingSystemRegions(sys_id, false);
			roi_info.genes = GlobalServiceProvider::database().processingSystemGenes(sys_id, false);
			config.roi = roi_info;
		}

		QString gsvar_xml_folder = Settings::path("gsvar_xml_folder", true);
		if (gsvar_xml_folder!="")
		{
			QString xml_file = gsvar_xml_folder + "/" + ps + "_tumor_only.xml" ;
			QByteArray temp_xml = Helper::tempFileName(".xml").toUtf8();
			worker.writeXML(temp_xml);
			ReportWorker::moveReport(temp_xml,xml_file);
		}
	}
	catch(Exception e)
	{
		QMessageBox::warning(this, "Could not create tumor only report", "Could not write tumor-only report. Error message: " + e.message());
	}

	QApplication::restoreOverrideCursor();

	//open report
	if (QMessageBox::question(this, "DNA tumor-only report", "report generated successfully!\nDo you want to open the report in your default RTF viewer?")==QMessageBox::Yes)
	{
		QDesktopServices::openUrl(QUrl::fromLocalFile(file_rep) );
	}
}


//transforms png data into list of tuples (png data in hex format, width, height)
QList<RtfPicture> pngsFromFiles(QStringList files)
{
	QList<RtfPicture> pic_list;
	foreach(const QString& path, files)
	{
		QImage pic;
		if (path.startsWith("http", Qt::CaseInsensitive))
		{
			QByteArray response = HttpHandler(true).get(path);
			if (!response.isEmpty()) pic.loadFromData(response);
		}
		else
		{
			pic = QImage(path);
		}
		if(pic.isNull()) continue;

		QByteArray png_data = "";
		QBuffer buffer(&png_data);
		buffer.open(QIODevice::WriteOnly);
		if (!pic.save(&buffer, "PNG")) continue;
		buffer.close();

		pic_list << RtfPicture(png_data.toHex(), pic.width(), pic.height());
	}

	return pic_list;
}

void MainWindow::generateReportSomaticRTF()
{
	if(!LoginManager::active()) return;

	NGSD db;
	QString ps_tumor = variants_.mainSampleName();
	QString ps_tumor_id = db.processedSampleId(ps_tumor);
	QString ps_normal = normalSampleName();
	QString ps_normal_id = db.processedSampleId(ps_normal);

	//Set data in somatic report settings
	somatic_report_settings_.report_config.setTargetRegionName(ui_.filters->targetRegion().name);

	somatic_report_settings_.report_config.setFilter((ui_.filters->filterName() != "[none]" ? ui_.filters->filterName() : "") ); //filter name -> goes to NGSD som. rep. conf.
	somatic_report_settings_.filters = ui_.filters->filters(); //filter cascase -> goes to report helper

	somatic_report_settings_.tumor_ps = ps_tumor;
	somatic_report_settings_.normal_ps = ps_normal;

	somatic_report_settings_.preferred_transcripts = GSvarHelper::preferredTranscripts();
	somatic_report_settings_.report_config.setEvaluationDate(QDate::currentDate());

	//load obo terms for filtering coding/splicing variants
	if (somatic_report_settings_.obo_terms_coding_splicing.size() == 0)
	{
		OntologyTermCollection obo_terms("://Resources/so-xp_3_1_0.obo",true);
		QList<QByteArray> ids;
		ids << obo_terms.childIDs("SO:0001580",true); //coding variants
		ids << obo_terms.childIDs("SO:0001568",true); //splicing variants
		foreach(const QByteArray& id, ids)
		{
			somatic_report_settings_.obo_terms_coding_splicing.add(obo_terms.getByID(id));
		}
	}

	somatic_report_settings_.target_region_filter = ui_.filters->targetRegion();
	if(!ui_.filters->targetRegion().isValid()) //use processing system data in case no filter is set
	{
		TargetRegionInfo generic_target;
		generic_target.regions = GlobalServiceProvider::database().processingSystemRegions(db.processingSystemIdFromProcessedSample(ps_tumor), false);
		generic_target.genes = db.genesToApproved(GlobalServiceProvider::database().processingSystemGenes(db.processingSystemIdFromProcessedSample(ps_tumor), false), true);
		generic_target.name = db.getProcessedSampleData(db.processedSampleId(ps_tumor)).processing_system;
		somatic_report_settings_.target_region_filter = generic_target;
	}

	if (db.getValues("SELECT value FROM processed_sample_qc AS psqc LEFT JOIN qc_terms as qc ON psqc.qc_terms_id = qc.id WHERE psqc.processed_sample_id=" + ps_tumor_id + " AND (qc.qcml_id ='QC:2000062' OR qc.qcml_id ='QC:2000063' OR qc.qcml_id ='QC:2000064') ").size() < 3)
	{
		QMessageBox::warning(this, "No HRD score found", "Warning:\nNo hrd score values found in the imported QC of tumor sample. HRD score set to 0.");
		somatic_report_settings_.report_config.setCnvLohCount(0);
		somatic_report_settings_.report_config.setCnvTaiCount(0);
		somatic_report_settings_.report_config.setCnvLstCount(0);
	}
	else
	{
		QString query = "SELECT value FROM processed_sample_qc AS psqc LEFT JOIN qc_terms as qc ON psqc.qc_terms_id = qc.id WHERE psqc.processed_sample_id=" + ps_tumor_id + " AND qc.qcml_id = :1";
		somatic_report_settings_.report_config.setCnvLohCount( db.getValue(query, false, "QC:2000062").toInt() );
		somatic_report_settings_.report_config.setCnvTaiCount( db.getValue(query, false, "QC:2000063").toInt() );
		somatic_report_settings_.report_config.setCnvLstCount( db.getValue(query, false, "QC:2000064").toInt() );
	}


	//Preselect report settings if not already exists to most common values
	if(db.somaticReportConfigId(ps_tumor_id, ps_normal_id) == -1)
	{
		somatic_report_settings_.report_config.setTumContentByMaxSNV(true);
		somatic_report_settings_.report_config.setTumContentByClonality(true);
		somatic_report_settings_.report_config.setTumContentByHistological(true);
		somatic_report_settings_.report_config.setMsiStatus(true);
		somatic_report_settings_.report_config.setCnvBurden(true);
	}

	//Parse genome ploidy from ClinCNV file
	FileLocation cnvFile = GlobalServiceProvider::fileLocationProvider().getAnalysisCnvFile();
	if (cnvFile.exists)
	{
		QStringList cnv_data = Helper::loadTextFile(cnvFile.filename, true, QChar::Null, true);

		for (const QString& line: cnv_data)
		{
			if (line.startsWith("##ploidy:"))
			{
				QStringList parts = line.split(':');
				somatic_report_settings_.report_config.setPloidy(parts[1].toDouble());
				break;
			}

			if (! line.startsWith("##"))
			{
				break;
			}
		}
	}

	//Get ICD10 diagnoses from NGSD
	QStringList tmp_icd10;
	QStringList tmp_phenotype;
	QStringList tmp_rna_ref_tissue;
	foreach(const auto& entry, db.getSampleDiseaseInfo(db.sampleId(ps_tumor)))
	{
		if(entry.type == "ICD10 code") tmp_icd10.append(entry.disease_info);
		if(entry.type == "clinical phenotype (free text)") tmp_phenotype.append(entry.disease_info);
		if(entry.type == "RNA reference tissue") tmp_rna_ref_tissue.append(entry.disease_info);
	}
	somatic_report_settings_.icd10 = tmp_icd10.join(", ");
	somatic_report_settings_.phenotype = tmp_phenotype.join(", ");

	SomaticReportDialog dlg(filename_, somatic_report_settings_, cnvs_, somatic_control_tissue_variants_, this); //widget for settings


	//Fill in RNA processed sample ids into somatic report dialog
	QSet<int> rna_ids =   db.relatedSamples(db.sampleId(ps_tumor).toInt(), "same sample", "RNA");
	if(!rna_ids.isEmpty())
	{
		dlg.enableChoiceRnaReportType(true);

		QStringList rna_names;
		foreach(int rna_id, rna_ids)
		{
			foreach(const auto& rna_ps_id, db.getValues("SELECT id FROM processed_sample WHERE sample_id=" + QString::number(rna_id)) )
			{
				rna_names << db.processedSampleName(rna_ps_id);
			}
		}
		dlg.setRNAids(rna_names);
	}

	// get all same samples
	int sample_id = db.sampleId(variants_.mainSampleName()).toInt();
	QSet<int> same_sample_ids = db.relatedSamples(sample_id, "same sample");
	same_sample_ids << sample_id; // add current sample id

	// get all related cfDNA
	QSet<int> cf_dna_sample_ids;
	foreach (int cur_sample_id, same_sample_ids)
	{
		cf_dna_sample_ids.unite(db.relatedSamples(cur_sample_id, "tumor-cfDNA"));
	}

	if (cf_dna_sample_ids.size() > 0)
	{
		dlg.enableChoicecfDnaReportType(true);
	}


	if(!dlg.exec())
	{
		return;
	}

	dlg.writeBackSettings();


	//store somatic report config in NGSD
	if(!dlg.skipNGSD())
	{
		db.setSomaticReportConfig(ps_tumor_id, ps_normal_id, somatic_report_settings_.report_config, variants_, cnvs_, somatic_control_tissue_variants_, Helper::userName());
	}

	QString destination_path; //path to rtf file
	if(dlg.getReportType() == SomaticReportDialog::report_type::DNA)
	{
		destination_path = last_report_path_ + "/" + ps_tumor + "_DNA_report_somatic_" + QDate::currentDate().toString("yyyyMMdd") + ".rtf";
	}
	else if (dlg.getReportType() == SomaticReportDialog::report_type::RNA)
	{
		destination_path = last_report_path_ + "/" + dlg.getRNAid() + "-" + ps_tumor + "_RNA_report_somatic_" + QDate::currentDate().toString("yyyyMMdd") + ".rtf";
	}
	else
	{
		destination_path = last_report_path_ + "/" + ps_tumor + "_cfDNA_report_somatic_" + QDate::currentDate().toString("yyyyMMdd") + ".rtf";
	}

	//get RTF file name
	QString file_rep = QFileDialog::getSaveFileName(this, "Store report file", destination_path, "RTF files (*.rtf);;All files(*.*)");
	if (file_rep=="") return;

	QApplication::setOverrideCursor(Qt::BusyCursor);

	if(dlg.getReportType() == SomaticReportDialog::report_type::DNA)
	{
		//generate somatic DNA report
		try
		{

			if(!SomaticReportHelper::checkGermlineSNVFile(somatic_control_tissue_variants_))
			{
				QApplication::restoreOverrideCursor();
				QMessageBox::warning(this, "Somatic report", "DNA report cannot be created because germline GSVar-file is invalid. Please check control tissue variant file.");
				return;
			}

			SomaticReportHelper report(GSvarHelper::build(), variants_, cnvs_, somatic_control_tissue_variants_, somatic_report_settings_);

			//Store XML file with the same somatic report configuration settings
			QTime timer;

			try
			{
				timer.start();
				QString tmp_xml = Helper::tempFileName(".xml");
				report.storeXML(tmp_xml);
				ReportWorker::moveReport(tmp_xml, Settings::path("gsvar_xml_folder") + "\\" + somatic_report_settings_.tumor_ps + "-" + somatic_report_settings_.normal_ps + ".xml");

				Log::perf("Generating somatic report XML took ", timer);
			}
			catch(Exception e)
			{
				QMessageBox::warning(this, "creation of XML file failed", e.message());
			}

			//Generate RTF
			timer.start();
			QByteArray temp_filename = Helper::tempFileName(".rtf").toUtf8();
			report.storeRtf(temp_filename);
			ReportWorker::moveReport(temp_filename, file_rep);
			Log::perf("Generating somatic report RTF took ", timer);

			//Generate files for QBIC upload
			timer.start();
			QString path = ps_tumor + "-" + ps_normal;
			if (GlobalServiceProvider::fileLocationProvider().isLocal()) path = Settings::string("qbic_data_path") + "/" + path;
			report.storeQbicData(path);
			Log::perf("Generating somatic report QBIC data took ", timer);

			QApplication::restoreOverrideCursor();
		}
		catch(Exception& error)
		{
			QApplication::restoreOverrideCursor();
			QMessageBox::warning(this, "Error while creating report", error.message());
			return;
		}
		catch(...)
		{
			QApplication::restoreOverrideCursor();
			QMessageBox::warning(this, "Error while creating report", "No error message!");
			return;
		}

		//open report
		if (QMessageBox::question(this, "DNA report", "DNA report generated successfully!\nDo you want to open the report in your default RTF viewer?")==QMessageBox::Yes)
		{
			QDesktopServices::openUrl(QUrl::fromLocalFile(file_rep) );
		}

		//reminder of MTB upload
		QStringList studies = db.getValues("SELECT s.name FROM study s, study_sample ss WHERE s.id=ss.study_id AND ss.processed_sample_id=" + ps_tumor_id);
		if (studies.contains("MTB"))
		{
			if (QMessageBox::question(this, "DNA report", "This sample is part of the study 'MTB'.\nDo you want to upload the data to MTB now?")==QMessageBox::Yes)
			{
				transferSomaticData();
			}
		}
	}
	else if (dlg.getReportType() == SomaticReportDialog::report_type::RNA)//RNA report
	{
		//Generate RTF
		try
		{
			QByteArray temp_filename = Helper::tempFileName(".rtf").toUtf8();

			SomaticRnaReportData rna_report_data = somatic_report_settings_;
			rna_report_data.rna_ps_name = dlg.getRNAid();
			rna_report_data.rna_fusion_file = GlobalServiceProvider::database().processedSamplePath(db.processedSampleId(dlg.getRNAid()), PathType::FUSIONS).filename;
			rna_report_data.rna_expression_file = GlobalServiceProvider::database().processedSamplePath(db.processedSampleId(dlg.getRNAid()), PathType::EXPRESSION).filename;
			rna_report_data.rna_bam_file = GlobalServiceProvider::database().processedSamplePath(db.processedSampleId(dlg.getRNAid()), PathType::BAM).filename;
			rna_report_data.ref_genome_fasta_file = Settings::string("reference_genome");

			try
			{
				QSharedPointer<VersatileFile> corr_file =  Helper::openVersatileFileForReading( GlobalServiceProvider::database().processedSamplePath( db.processedSampleId(dlg.getRNAid()), PathType::EXPRESSION_CORR ).filename );
				rna_report_data.expression_correlation = Helper::toDouble(corr_file->readAll());
			}
			catch(Exception)
			{
				rna_report_data.expression_correlation = std::numeric_limits<double>::quiet_NaN();
			}

			try
			{
				TSVFileStream cohort_file( GlobalServiceProvider::database().processedSamplePath( db.processedSampleId(dlg.getRNAid()), PathType::EXPRESSION_COHORT ).filename );
				rna_report_data.cohort_size = cohort_file.header().count()-1;
			}
			catch(Exception)
			{
			}

			rna_report_data.rna_qcml_data = db.getQCData(db.processedSampleId(dlg.getRNAid()));

			//Add data from fusion pics
			try
			{
				rna_report_data.fusion_pics = pngsFromFiles(GlobalServiceProvider::database().getRnaFusionPics(dlg.getRNAid()));
			}
			catch(Exception) //Nothing to do here
			{
			}
			//Add data from expression plots
			try
			{
				rna_report_data.expression_plots = pngsFromFiles(GlobalServiceProvider::database().getRnaExpressionPlots(dlg.getRNAid()));
			}
			catch(Exception)
			{
			}

			//Look in tumor sample for HPA reference tissue
			foreach(const auto& entry, db.getSampleDiseaseInfo(db.sampleId(dlg.getRNAid())) )
			{
				if(entry.type == "RNA reference tissue") tmp_rna_ref_tissue.append(entry.disease_info);
			}
			tmp_rna_ref_tissue.removeDuplicates();
			rna_report_data.rna_hpa_ref_tissue = tmp_rna_ref_tissue.join(", ");

			SomaticRnaReport rna_report(variants_, cnvs_, rna_report_data);

			rna_report.writeRtf(temp_filename);
			ReportWorker::moveReport(temp_filename, file_rep);
			QApplication::restoreOverrideCursor();
		}
		catch(Exception& error)
		{
			QApplication::restoreOverrideCursor();
			QMessageBox::warning(this, "Error while creating somatic RNA report.", error.message());
			return;
		}
		catch(...)
		{
			QApplication::restoreOverrideCursor();
			QMessageBox::warning(this, "Error while creating somatic RNA report.", "No error message!");
			return;
		}

		if (QMessageBox::question(this, "RNA report", "RNA report generated successfully!\nDo you want to open the report in your default RTF viewer?")==QMessageBox::Yes)
		{
			QDesktopServices::openUrl(QUrl::fromLocalFile(file_rep));
		}
	}
	else if (dlg.getReportType() == SomaticReportDialog::report_type::cfDNA)
	{
		try
		{
			QStringList errors;
			CfdnaDiseaseCourseTable table = GSvarHelper::cfdnaTable(ps_tumor, errors, true);
			SomaticcfDNAReportData data(somatic_report_settings_, table);
			SomaticcfDnaReport report(data);

			QByteArray temp_filename = Helper::tempFileName(".rtf").toUtf8();
			report.writeRtf(temp_filename);
			ReportWorker::moveReport(temp_filename, file_rep);
			QApplication::restoreOverrideCursor();

			if (QMessageBox::question(this, "cfDNA report", "cfDNA report generated successfully!\nDo you want to open the report in your default RTF viewer?")==QMessageBox::Yes)
			{
				QDesktopServices::openUrl(QUrl::fromLocalFile(file_rep));
			}
		}
		catch (Exception& error)
		{
			QApplication::restoreOverrideCursor();
			QMessageBox::warning(this, "Error while gathering data for somatic cfDNA report.", error.message());
			return;
		}
	}
	else
	{
		QApplication::restoreOverrideCursor();
		QMessageBox::warning(this, "Unknown somatic report type!", "Unknown somatic report type! This should not happen please inform the bioinformatic team.");
		return;
	}
}

void MainWindow::generateReportGermline()
{
	//check that sample is in NGSD
	NGSD db;
	QString ps_name = germlineReportSample();
	QString sample_id = db.sampleId(ps_name, false);
	QString processed_sample_id = db.processedSampleId(ps_name, false);
	if (sample_id.isEmpty() || processed_sample_id.isEmpty())
	{
		GUIHelper::showMessage("Error", "Sample not found in the NGSD.\nCannot generate a report for samples that are not in the NGSD!");
		return;
	}

	//check if there are unclosed gaps
	QStringList unclosed_gap_ids = db.getValues("SELECT id FROM gaps WHERE processed_sample_id='" + processed_sample_id + "' AND (status='to close' OR status='in progress')");
	if (unclosed_gap_ids.count()>0 && QMessageBox::question(this, "Not all gaps closed", "There are gaps for this sample, which still have to be closed!\nDo you want to continue?")==QMessageBox::No)
	{
		return;
	}

	//show report dialog
	ReportDialog dialog(ps_name, report_settings_, variants_, cnvs_, svs_, ui_.filters->targetRegion(), this);
	if (!dialog.exec()) return;

	//set report type
	report_settings_.report_type = dialog.type();

	//get export file name
	QString trio_suffix = (variants_.type() == GERMLINE_TRIO ? "trio_" : "");
	QString type_suffix = dialog.type();
	if (type_suffix!="all") type_suffix = type_suffix.replace(" ", "_") + "s";
	QString roi_name = ui_.filters->targetRegion().name;
	if (roi_name!="") //remove date and prefix with '_'
	{
		roi_name.remove(QRegExp("_[0-9]{4}_[0-9]{2}_[0-9]{2}"));
		roi_name = "_" + roi_name;
	}
	QString file_rep = QFileDialog::getSaveFileName(this, "Export report file", last_report_path_ + "/" + ps_name + roi_name + "_report_" + trio_suffix + type_suffix + "_" + QDate::currentDate().toString("yyyyMMdd") + ".html", "HTML files (*.html);;All files(*.*)");
	if (file_rep=="") return;
	last_report_path_ = QFileInfo(file_rep).absolutePath();

	//prepare report generation data
	PrsTable prs_table;
	FileLocationList prs_files = GlobalServiceProvider::fileLocationProvider().getPrsFiles(false).filterById(ps_name);
	if (prs_files.count()==1) prs_table.load(prs_files[0].filename);

	GermlineReportGeneratorData data(GSvarHelper::build(), ps_name, variants_, cnvs_, svs_, prs_table, report_settings_, ui_.filters->filters(), GSvarHelper::preferredTranscripts(), GlobalServiceProvider::statistics());
	data.processing_system_roi = GlobalServiceProvider::database().processingSystemRegions(db.processingSystemIdFromProcessedSample(ps_name), false);
	data.ps_bam = GlobalServiceProvider::database().processedSamplePath(processed_sample_id, PathType::BAM).filename;
	data.ps_lowcov = GlobalServiceProvider::database().processedSamplePath(processed_sample_id, PathType::LOWCOV_BED).filename;
	if (ui_.filters->targetRegion().isValid())
	{
		data.roi = ui_.filters->targetRegion();
		data.roi.genes = db.genesToApproved(data.roi.genes, true);
	}

	//show busy dialog
	busy_dialog_ = new BusyDialog("Report", this);
	busy_dialog_->init("Generating report", false);

	//start worker in new thread
	ReportWorker* worker = new ReportWorker(data, file_rep);
	connect(worker, SIGNAL(finished(bool)), this, SLOT(reportGenerationFinished(bool)));
	worker->start();
}

void MainWindow::reportGenerationFinished(bool success)
{
	delete busy_dialog_;

	//show result info box
	ReportWorker* worker = qobject_cast<ReportWorker*>(sender());
	if (success)
	{
		if (QMessageBox::question(this, "Report", "Report generated successfully!\nDo you want to open the report in your browser?")==QMessageBox::Yes)
		{
			QDesktopServices::openUrl(QUrl::fromLocalFile(worker->getReportFile()));
		}
	}
	else
	{
		QMessageBox::warning(this, "Error", "Report generation failed:\n" + worker->errorMessage());
	}
	//clean
	worker->deleteLater();
}

void MainWindow::openProcessedSampleTabsCurrentAnalysis()
{
	if (filename_=="") return;

	SampleHeaderInfo infos = variants_.getSampleHeader();
	foreach(const SampleInfo& info, infos)
	{
		openProcessedSampleTab(info.id);
	}
}

void MainWindow::on_actionOpenProcessedSampleTabByName_triggered()
{
	ProcessedSampleSelector dlg(this, false);
	if (!dlg.exec()) return;

	QString ps_name = dlg.processedSampleName();
	if (ps_name.isEmpty()) return;

	openProcessedSampleTab(ps_name);
}

void MainWindow::on_actionOpenSequencingRunTabByName_triggered()
{
	//create
	DBSelector* selector = new DBSelector(this);
	NGSD db;
	selector->fill(db.createTable("sequencing_run", "SELECT id, name FROM sequencing_run"));

	//show
	auto dlg = GUIHelper::createDialog(selector, "Select sequencing run", "run:", true);
	if (dlg->exec()==QDialog::Rejected) return ;

	//handle invalid name
	if (selector->getId()=="") return;

	openRunTab(selector->text());
}

QString MainWindow::selectGene()
{
	//create
	DBSelector* selector = new DBSelector(QApplication::activeWindow());
	NGSD db;
	selector->fill(db.createTable("gene", "SELECT id, symbol FROM gene"));

	//show
	auto dlg = GUIHelper::createDialog(selector, "Select gene", "symbol (or transcript name):", true);
	if (dlg->exec()==QDialog::Rejected) return "";

	//handle invalid gene name > check if it is a transcript name
	if (selector->getId()=="")
	{
		int gene_id = db.geneIdOfTranscript(selector->text().toUtf8(), false, GSvarHelper::build());
		if (gene_id!=-1)
		{
			return db.geneSymbol(gene_id);
		}
	}

	return selector->text();
}

QString MainWindow::selectProcessedSample()
{
	//determine processed sample names
	QStringList ps_list;
	foreach(const SampleInfo& info, variants_.getSampleHeader())
	{
		ps_list << info.id.trimmed();
	}

	//no samples => error
	if (ps_list.isEmpty())
	{
		THROW(ProgrammingException, "selectProcessedSample() cannot be used if there is no variant list loaded!");
	}

	//one sample => auto-select
	if (ps_list.count()==1)
	{
	   return ps_list[0];
	}

	//several affected => let user select
	bool ok = false;
	QString selected = QInputDialog::getItem(this, "Select processed sample", "processed sample:", ps_list, 0, false, &ok);
	if (ok) return selected;

	return "";
}

const TargetRegionInfo& MainWindow::targetRegion()
{
	return ui_.filters->targetRegion();
}

void MainWindow::importBatch(QString title, QString text, QString table, QStringList fields)
{
	//show dialog
	QTextEdit* edit = new QTextEdit();
	edit->setAcceptRichText(false);
	auto dlg = GUIHelper::createDialog(edit, title, text, true);
	if (dlg->exec()!=QDialog::Accepted) return;

	// load input text as table
	QList<QStringList> table_content;
	QStringList lines = edit->toPlainText().split("\n");
	foreach(const QString& line, lines)
	{
		// skip empty lines
		if (line.trimmed().isEmpty()) continue;

		table_content.append(line.split("\t"));
	}

	NGSD db;

	//special handling of processed sample: add 'process_id' and get the list of all cfDNA processing systems
	QStringList cfdna_processing_systems;
	if (table=="processed_sample")
	{
		fields.append("process_id");

		// get all cfDNA processing system
		cfdna_processing_systems = db.getValues("SELECT name_manufacturer FROM processing_system WHERE type='cfDNA (patient-specific)'");
	}

	// special handling of sample: extract last column containing the corresponding tumor sample id for a cfDNA sample
	QList<QPair<QString,QString>> cfdna_tumor_relation;
	if (table=="sample")
	{
		int name_idx = fields.indexOf("name");
		for (int i = 0; i < table_content.size(); ++i)
		{
			QStringList& row = table_content[i];
			if (row.length() == (fields.length() + 1))
			{
				//store tumor-cfDNA relation
				QString tumor_sample = row.last();
				QString cfdna_sample = row.at(name_idx).trimmed();
				// only import if not empty
				if (!tumor_sample.isEmpty()) cfdna_tumor_relation.append(QPair<QString,QString>(tumor_sample, cfdna_sample));

				// remove last element
				row.removeLast();
			}
		}
	}

	//special handling of sample_relations: add user
	if (table=="sample_relations")
	{
		fields.append("user_id");
	}

	//special handling of sample_disease_info: add type and user
	if (table=="sample_disease_info")
	{
		fields.append("type");
		fields.append("user_id");
	}

	//prepare query
	QString query_str = "INSERT INTO " + table + " (" + fields.join(", ") + ") VALUES (";
	for(int i=0; i<fields.count(); ++i)
	{
		if (i!=0) query_str += ", ";
		query_str += ":" + QString::number(i);
	}
	query_str += ")";

	SqlQuery q_insert = db.getQuery();
	q_insert.prepare(query_str);

	//check and insert
	QString last_processed_line;
	try
	{
		db.transaction();

		int imported = 0;
		QStringList duplicates;
		QStringList missing_cfDNA_relation;
		for (int i = 0; i < table_content.size(); ++i)
		{
			QStringList& row = table_content[i];
			last_processed_line = row.join("\t");

			//special handling of processed sample: add 'process_id'
			if (table=="processed_sample")
			{
				QString sample_name = row[0].trimmed();
				QString sample_id = db.sampleId(sample_name);

				QString next_id = db.nextProcessingId(sample_id);
				if (next_id.toInt() > 99)
				{
					THROW(ArgumentException, "Error: For sample " + sample_name + " already exist 99 processed samples.\nCannot create more processed samples for " + sample_name + ".\n\nPlease import a new sample ex. NA12878 -> NA12878x2");
				}
				row.append(next_id);
			}

			//special handling of sample duplicates
			if (table=="sample")
			{
				QString sample_name = row[0].trimmed();
				QString sample_id = db.sampleId(sample_name, false);
				if (sample_id!="")
				{
					duplicates << sample_name;
					continue;
				}
			}

			//special handling of sample_relations: add user
			if (table=="sample_relations")
			{
				row.append(LoginManager::userName());
			}

			//special handling of sample_disease_info:
			if (table=="sample_disease_info")
			{
				//add type and user
				row.append("HPO term id");
				row.append(LoginManager::userName());

				//skip duplicates
				QString sample_name = row[0].trimmed();
				QString hpo_id = row[1].trimmed();
				QString sample_id = db.sampleId(sample_name, false);
				QList<int> duplicate_ids = db.getValuesInt("SELECT id FROM sample_disease_info WHERE sample_id='"+sample_id+"' AND disease_info='"+hpo_id+"' AND type='HPO term id'");
				if (!duplicate_ids.isEmpty())
				{
					duplicates << sample_name+"/"+hpo_id;
					continue;
				}
			}

			//special handling of study-sample relation
			if (table=="study_sample")
			{
				//skip duplicates
				QString study = row[0].trimmed();
				QString study_id = db.getValue("SELECT id FROM study WHERE name=:0", true, study).toString();
				QString ps = row[1].trimmed();
				QString ps_id = db.processedSampleId(ps, false);
				QList<int> duplicate_ids = db.getValuesInt("SELECT id FROM study_sample ss WHERE processed_sample_id='"+ps_id+"' AND study_id='"+study_id+"'");
				if (!duplicate_ids.isEmpty())
				{
					duplicates << study+"/"+ps;
					continue;
				}
			}

			//check tab-separated parts count
			if (row.count()!=fields.count()) THROW(ArgumentException, "Error: line with more/less than " + QString::number(fields.count()) + " tab-separated parts.");
			//check and bind
			for(int i=0; i<fields.count(); ++i)
			{
				QString field = fields[i];
				const TableFieldInfo& field_info = db.tableInfo(table).fieldInfo(field);

				QString value = row[i].trimmed();

				//check for cfDNA samples if a corresponding tumor sample exists
				if (table=="processed_sample" && field=="processing_system_id" && cfdna_processing_systems.contains(value))
				{
					// get tumor relation
					int sample_id = db.sampleId(row[0].trimmed()).toInt();
					if (db.relatedSamples(sample_id, "tumor-cfDNA").size() < 1)
					{
						//No corresponding tumor found!
						missing_cfDNA_relation.append(row[0].trimmed());
					}
				}

				//special handling of sample_disease_info
				if (table=="sample_disease_info" && field=="disease_info")
				{
					//check HPO term id is valid
					db.getValue("SELECT id FROM hpo_term WHERE hpo_id=:0", false, value);
				}

				//FK: name to id
				if (field_info.type==TableFieldInfo::FK && !value.isEmpty())
				{
					QString name_field = field_info.fk_name_sql;
					if (name_field.startsWith("CONCAT(name")) //some FK-fields show additional information after the name > use only the name
					{
						name_field = name_field.left(name_field.indexOf(','));
						name_field = name_field.mid(7);
					}
					value = db.getValue("SELECT " + field_info.fk_field + " FROM " + field_info.fk_table + " WHERE " + name_field + "=:0", false, value).toString();
				}

				//accept German dates as well
				if (field_info.type==TableFieldInfo::DATE && !value.isEmpty())
				{
					QDate date_german = QDate::fromString(value, "dd.MM.yyyy");
					if (date_german.isValid())
					{
						value = date_german.toString(Qt::ISODate);
					}
				}

				//decimal point for numbers
				if (field_info.type==TableFieldInfo::FLOAT && value.contains(','))
				{
					value = value.replace(',', '.');
				}

				//check errors
				QStringList errors = db.checkValue(table, field, value, true);
				if (errors.count()>0)
				{
					THROW(ArgumentException, "Error: Invalid value '" + value + "' for field '" + field + "':\n" + errors.join("\n"));
				}

				q_insert.bindValue(i, value.isEmpty() && field_info.is_nullable ? QVariant() : value);
			}

			//insert
			q_insert.exec();
			++imported;
		}

		//ask user if duplicates should be skipped
		if (duplicates.count()>0)
		{
			int button = QMessageBox::question(this, title, QString::number(duplicates.count()) +" entries are already present in the NGSD:\n" + duplicates.join(" ") + "\n\n Do you want to skip these duplicates and continue?", QMessageBox::Ok, QMessageBox::Abort);
			if (button==QMessageBox::Abort)
			{
				db.rollback();
				return;
			}
		}

		//import tumor-cfDNA relations
		int imported_relations = 0;
		for (int i = 0; i < cfdna_tumor_relation.size(); ++i)
		{
			const QPair<QString, QString>& relation = cfdna_tumor_relation.at(i);
			// check for valid input
			if (db.getSampleData(db.sampleId(relation.first)).is_tumor)
			{
				// add relation
				db.addSampleRelation(SampleRelation{relation.first.toUtf8(), "tumor-cfDNA", relation.second.toUtf8()}, true);
				imported_relations++;
			}
			else
			{
				THROW(DatabaseException, "Sample " + relation.first + " is not a tumor! Can't import relation.");
			}
		}

		// abort if cfDNA-tumor relations are missing
		if (missing_cfDNA_relation.size() > 0)
		{
			THROW(ArgumentException, "The NGSD does not contain a cfDNA-tumor relation for the following cfDNA samples:\n" + missing_cfDNA_relation.join(", ") + "\nAborting import.");
		}

		db.commit();

		QString message = "Imported " + QString::number(imported) + " table rows.";
		if (duplicates.count()>0) message += "\n\nSkipped " + QString::number(duplicates.count()) + " table rows.";
		if (imported_relations > 0) message += "\n\nImported " + QString::number(imported_relations) + " tumor-cfDNA relations.";
		QMessageBox::information(this, title,  message);
	}
	catch (Exception& e)
	{
		db.rollback();
		QMessageBox::warning(this, title + " - failed", "Import failed - no data was imported!\n\nLine:\n"+ last_processed_line.trimmed() + "\n\nError message:\n" + e.message());
	}
}

const VariantList&MainWindow::getSmallVariantList()
{
	return variants_;
}

const CnvList&MainWindow::getCnvList()
{
	return cnvs_;
}

const BedpeFile&MainWindow::getSvList()
{
	return svs_;
}

void MainWindow::on_actionOpenGeneTabByName_triggered()
{
	QString symbol = selectGene();
	if (symbol=="") return;

	openGeneTab(symbol);
}


void MainWindow::on_actionOpenVariantTab_triggered()
{
	try
	{
		VariantOpenDialog dlg(this);

		if (dlg.exec()!=QDialog::Accepted) return;

		openVariantTab(dlg.variant());
	}
	catch(Exception& e)
	{
		QMessageBox::warning(this, "Invalid variant text", e.message());
	}
}

void MainWindow::on_actionOpenProcessingSystemTab_triggered()
{
	//create
	DBSelector* selector = new DBSelector(this);
	NGSD db;
	selector->fill(db.createTable("processing_system", "SELECT id, CONCAT(name_manufacturer, ' (', name_short, ')') FROM processing_system"));

	//show
	auto dlg = GUIHelper::createDialog(selector, "Select processing system", "name:", true);
	if (dlg->exec()==QDialog::Rejected) return;

	//handle invalid name
	if (selector->getId()=="") return;

	openProcessingSystemTab(db.getValue("SELECT name_short FROM processing_system WHERE id=" + selector->getId()).toString());
}

void MainWindow::on_actionOpenProjectTab_triggered()
{
	//create
	DBSelector* selector = new DBSelector(this);
	NGSD db;
	selector->fill(db.createTable("project", "SELECT id, name FROM project"));

	//show
	auto dlg = GUIHelper::createDialog(selector, "Select project", "project:", true);
	if (dlg->exec()==QDialog::Rejected) return ;

	//handle invalid name
	if (selector->getId()=="") return;

	openProjectTab(selector->text());
}

void MainWindow::on_actionStatistics_triggered()
{
	try
	{
		LoginManager::checkRoleIn(QStringList{"admin", "user"});
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Permissions error", e.message());
		return;
	}

	QApplication::setOverrideCursor(Qt::BusyCursor);

	NGSD db;
	TsvFile table;
	bool human_only = true;

	//table header
	table.addHeader("month");
	QStringList sys_types = QStringList() << "WGS" << "WES" << "Panel" << "RNA" << "lrGS";
	QStringList pro_types = QStringList() << "diagnostic" << "research";
	foreach(QString pro, pro_types)
	{
		foreach(QString sys, sys_types)
		{
			table.addHeader(sys + " " + pro);
		}
	}

	//table rows
	QSet<QString> comments;
	QDate start = QDate::currentDate();
	start = start.addDays(1-start.day());
	QDate end = start.addMonths(1);
	while(start.year()>=2015)
	{
		QVector<int> counts(table.headers().count(), 0);

		//select runs of current month
		SqlQuery q_run_ids = db.getQuery();
		q_run_ids.exec("SELECT id FROM sequencing_run WHERE end_date >='" + start.toString(Qt::ISODate) + "' AND end_date < '" + end.toString(Qt::ISODate) + "' AND status!='analysis_not_possible' AND status!='run_aborted' AND status!='n/a' AND quality!='bad'");
		while(q_run_ids.next())
		{
			//select samples
			SqlQuery q_sample_data = db.getQuery();
			q_sample_data.exec("SELECT sys.type, p.type FROM sample s, processed_sample ps, processing_system sys, project p WHERE " + QString(human_only ? " s.species_id=(SELECT id FROM species WHERE name='human') AND " : "") + " ps.processing_system_id=sys.id AND ps.sample_id=s.id AND ps.project_id=p.id AND ps.sequencing_run_id=" + q_run_ids.value(0).toString() + " AND ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples)");

			//count
			while(q_sample_data.next())
			{
				QString sys_type = q_sample_data.value(0).toString();
				if (sys_type.contains("Panel")) sys_type = "Panel";
				if (!sys_types.contains(sys_type))
				{
					comments << "##Skipped processing system type '" + sys_type + "'";
					continue;
				}
				QString pro_type = q_sample_data.value(1).toString();
				if (!pro_types.contains(pro_type))
				{
					comments << "##Skipped project type '" + pro_type + "'";
					continue;
				}

				int index = table.headers().indexOf(sys_type + " " + pro_type);
				++counts[index];
			}
		}

		//create row
		QStringList row;
		for(int i=0; i<counts.count(); ++i)
		{
			if (i==0)
			{
				row << start.toString("yyyy/MM");
			}
			else
			{
				row << QString::number(counts[i]);
			}
		}
		table.addRow(row);

		//next month
		start = start.addMonths(-1);
		end = end.addMonths(-1);
	}

	//comments
	foreach(QString comment, comments)
	{
		table.addComment(comment);
	}

	QApplication::restoreOverrideCursor();

	//show dialog
	TsvTableWidget* widget = new TsvTableWidget(table);
	widget->setMinimumWidth(850);
	auto dlg = GUIHelper::createDialog(widget, "Statistics", "Sequencing statistics grouped by month (human)");
	dlg->exec();
}

void MainWindow::on_actionDevice_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("device");
	auto dlg = GUIHelper::createDialog(widget, "Device administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionGenome_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("genome");
	auto dlg = GUIHelper::createDialog(widget, "Genome administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionMID_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("mid");
	auto dlg = GUIHelper::createDialog(widget, "MID administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionProcessingSystem_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("processing_system");
	auto dlg = GUIHelper::createDialog(widget, "Processing system administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionProject_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("project");
	auto dlg = GUIHelper::createDialog(widget, "Project administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionSample_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("sample");
	auto dlg = GUIHelper::createDialog(widget, "Sample administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionSampleGroup_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("sample_group");
	auto dlg = GUIHelper::createDialog(widget, "Sample group administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionSender_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("sender");
	auto dlg = GUIHelper::createDialog(widget, "Sender administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionSpecies_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("species");
	auto dlg = GUIHelper::createDialog(widget, "Species administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionUsers_triggered()
{
	try
	{
		LoginManager::checkRoleIn(QStringList{"admin"});
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Permissions error", e.message());
		return;
	}

	//show user table
	DBTableAdministration* widget = new DBTableAdministration("user");
	auto dlg = GUIHelper::createDialog(widget, "User administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionExportTestData_triggered()
{
	NGSD db;
	QMap<QString, QSet<int>> sql_history;
	QTime timer;
	QStringList base_tables = {
		"user",
		"device",
		"disease_term",
		"disease_gene",
		"gene",
		"gene_alias",
		"gene_transcript",
		"gene_exon",
		"gene_pseudogene_relation",
		"geneinfo_germline",
		"genome",
		"hpo_term",
		"hpo_parent",
		"hpo_genes",
		"mid",
		"omim_gene",
		"omim_phenotype",
		"omim_preferred_phenotype",
		"preferred_transcripts",
		"processing_system",
		"project",
		"qc_terms",
		"sender",
		"sequencing_run",
		"somatic_pathway",
		"somatic_pathway_gene",
		"somatic_gene_role",
		"runqc_read",
		"runqc_lane",
		"species"
	};

	try
	{
		LoginManager::checkRoleIn(QStringList{"admin", "user"});

		//get and check processed sample list
		QStringList ps_list;
		QString ps_text = QInputDialog::getMultiLineText(this, "Test data export", "List the processed samples (one per line):");
		foreach(const QString& ps, ps_text.split("\n"))
		{
			QString ps_id = db.processedSampleId(ps);
			QString project_type = db.getProcessedSampleData(ps_id).project_type;
			if (project_type!="test")
			{
				THROW(ArgumentException, "Processes sample '" + ps + "' has project type '" +project_type + "', but only samples from type 'test' can be exported!");
			}
			ps_list << ps;
		}
		if (ps_list.count() == 0) return;

		//get and open output file
		QString file_name = QFileDialog::getSaveFileName(this, "Export database tables", QDir::homePath()+QDir::separator()+"db_data_"+QDateTime::currentDateTime().toString("dd_MM_yyyy")+".sql", "SQL (*.sql);;All files (*.*)");
		if (file_name.isEmpty()) return;

		QSharedPointer<QFile> file = Helper::openFileForWriting(file_name, false);
		QTextStream output_stream(file.data());

		QApplication::setOverrideCursor(Qt::BusyCursor);

		timer.start();
		for (int i = 0; i < base_tables.count(); i++)
		{
			ui_.statusBar->showMessage("Exporting table \"" + base_tables[i] + "\"");
			QApplication::processEvents();
			db.exportTable(base_tables[i], output_stream);
		}
		Log::perf("Exporting base tables took ", timer);

		timer.start();
		foreach(const QString& ps, ps_list)
		{
			ui_.statusBar->showMessage("Exporting data of " + ps);
			QApplication::processEvents();

			QString s_id = db.sampleId(ps);
			QString ps_id = db.processedSampleId(ps);
			db.exportTable("sample", output_stream, "id='"+s_id+"'", &sql_history);
			db.exportTable("sample_disease_info", output_stream, "sample_id='"+s_id+"'", &sql_history);
			db.exportTable("processed_sample", output_stream, "id='"+ps_id+"'", &sql_history);
			db.exportTable("processed_sample_qc", output_stream, "processed_sample_id='"+ps_id+"'", &sql_history);

			QStringList variant_id_list = db.getValues("SELECT variant_id FROM detected_variant WHERE processed_sample_id='"+ps_id+"'");
			db.exportTable("variant", output_stream, "id IN ("+variant_id_list.join(", ")+")", &sql_history);
			db.exportTable("detected_variant", output_stream, "processed_sample_id='"+ps_id+"'", &sql_history);

			QString ps_cnv_id = db.getValue("SELECT id FROM cnv_callset WHERE processed_sample_id=:0", true, ps_id).toString();
			if (!ps_cnv_id.isEmpty())
			{
				db.exportTable("cnv_callset", output_stream, "id="+ps_cnv_id, &sql_history);
				db.exportTable("cnv", output_stream, "cnv_callset_id="+ps_cnv_id, &sql_history);
			}

			QString sv_callset_id = db.getValue("SELECT id FROM sv_callset WHERE processed_sample_id='"+ps_id+"'", true).toString();
			if (!sv_callset_id.isEmpty())
			{
				db.exportTable("sv_callset", output_stream, "id="+sv_callset_id, &sql_history);
				db.exportTable("sv_deletion", output_stream, "sv_callset_id="+sv_callset_id, &sql_history);
				db.exportTable("sv_duplication", output_stream, "sv_callset_id="+sv_callset_id, &sql_history);
				db.exportTable("sv_insertion", output_stream, "sv_callset_id="+sv_callset_id, &sql_history);
				db.exportTable("sv_inversion", output_stream, "sv_callset_id="+sv_callset_id, &sql_history);
				db.exportTable("sv_translocation", output_stream, "sv_callset_id="+sv_callset_id, &sql_history);
			}
		}
		Log::perf("Exporting processed sample data took ", timer);

		QApplication::restoreOverrideCursor();

		QMessageBox::information(this, "Test data export", "Exported test data to " + file_name);
	}
	catch (Exception& e)
	{
		GUIHelper::showException(this, e, "NGSD export error");
	}
}

void MainWindow::on_actionImportSequencingRuns_triggered()
{
	importBatch("Import sequencing runs",
				"Batch import of sequencing runs. Must contain the following tab-separated fields:<br><b>name</b>, flowcell ID, <b>flowcell type</b>, start date, end date, <b>device</b>, <b>recipe</b>, pool_molarity, <b>pool quantification method</b>, comment",
				"sequencing_run",
				QStringList() << "name" << "fcid" << "flowcell_type" << "start_date" << "end_date" << "device_id" << "recipe" << "pool_molarity" << "pool_quantification_method" << "comment"
				);
}

void MainWindow::on_actionImportTestData_triggered()
{
	NGSD db;

	try
	{
		//check role
		LoginManager::checkRoleIn(QStringList{"admin", "user"});

		//prevent overriding the production database
		if (db.isProductionDb())
		{
			THROW(DatabaseException, "Cannot import test data into a production database!");
		}

		//get input file
		QString file_name = QFileDialog::getOpenFileName(this, "Import SQL data", QDir::homePath(), "SQL (*.sql);;All files (*.*)");
		if (file_name.isEmpty()) return;

		//import
		QApplication::setOverrideCursor(Qt::BusyCursor);
		db.removeInitData();

		QSharedPointer<QFile> file = Helper::openFileForReading(file_name, false);
		while(!file->atEnd())
		{
			QString line = file->readLine().trimmed();
			if (line.isEmpty()) continue;

			//comments > show in status bar
			if (line.startsWith("--"))
			{
				line = line.replace("--", "").trimmed();
				if (!line.isEmpty())
				{
					ui_.statusBar->showMessage("Importing \"" + line + "\"");
					QApplication::processEvents();
				}
				continue;
			}

			//import line
			db.getQuery().exec(line);
		}

		QApplication::restoreOverrideCursor();

		QMessageBox::information(this, "Test data import", "Import is complete");
	}
	catch (Exception& e)
	{
		GUIHelper::showException(this, e, "NGSD import error");
	}
}

void MainWindow::on_actionImportMids_triggered()
{
	importBatch("Import MIDs",
				"Batch import of MIDs. Please enter MIDs as tab-delimited text.<br>Example:<br><br>illumina 1 → CGTGAT<br>illumina 2 → AGATA<br>illumina 3 → GTCATG",
				 "mid",
				QStringList() << "name" << "sequence"
				);
}

void MainWindow::on_actionImportStudy_triggered()
{
	importBatch("Import study",
				"Batch import of stamples to studies. Please enter study, processed sample and study-specific name of sample (can be empty).<br>Example:<br><br>SomeStudy → NA12345_01 → NameOfSampleInStudy",
				 "study_sample",
				QStringList() << "study_id" << "processed_sample_id" << "study_sample_idendifier"
				);
}
void MainWindow::on_actionImportSamples_triggered()
{
	importBatch("Import samples",
				"Batch import of samples. Must contain the following tab-separated fields:<br><b>name</b>, name external, <b>sender</b>, received, received by, <b>sample type</b>, <b>tumor</b>, <b>ffpe</b>, <b>species</b>, concentration [ng/ul], volume, 260/280, 260/230, RIN/DIN, <b>gender</b>, <b>quality</b>, comment, disease group, disease status, tissue. <br> (For cfDNA Samples a additional column which defines the corresponding tumor sample can be given.) ",
				"sample",
				QStringList() << "name" << "name_external" << "sender_id" << "received" << "receiver_id" << "sample_type" << "tumor" << "ffpe" << "species_id" << "concentration" << "volume" << "od_260_280" << "od_260_230" << "integrity_number" << "gender" << "quality" << "comment" << "disease_group" << "disease_status" << "tissue"
				);
}

void MainWindow::on_actionImportProcessedSamples_triggered()
{
	importBatch("Import processed samples",
				"Batch import of processed samples. Must contain the following tab-separated fields:<br><b>sample</b>, <b>project</b>, <b>run name</b>, <b>lane</b>, mid1 name, mid2 name, operator, <b>processing system</b>, processing input [ng], molarity [nM], comment, normal processed sample, <b>processing modus</b>, batch number",
				"processed_sample",
				QStringList() << "sample_id" << "project_id" << "sequencing_run_id" << "lane" << "mid1_i7" << "mid2_i5" << "operator_id" << "processing_system_id" << "processing_input" << "molarity" << "comment" << "normal_id" << "processing_modus" << "batch_number"
				);
}

void MainWindow::on_actionImportSampleRelations_triggered()
{
	importBatch("Import sample relations",
				"Batch import of sample relations. Must contain the following tab-separated fields:<br><b>sample1</b>, <b>relation</b>, <b>sample2</b>",
				"sample_relations",
				QStringList() << "sample1_id" << "relation" << "sample2_id"
				);
}

void MainWindow::on_actionImportSampleHpoTerms_triggered()
{
	importBatch("Import sample HPO terms",
				"Batch import of sample HPO terms. Must contain the following tab-separated fields:<br><b>sample1</b>, <b>HPO term id e.g. 'HP:0003002'</b>",
				"sample_disease_info",
				QStringList() << "sample_id" << "disease_info"
				);
}

void MainWindow::on_actionImportCfDNAPanels_triggered()
{
	CfDNAPanelBatchImport* widget = new CfDNAPanelBatchImport();
//	auto dlg = GUIHelper::createDialog(widget, "Import cfDNA panels");
//	addModelessDialog(dlg);
	widget->exec();
}

void MainWindow::on_actionMidClashDetection_triggered()
{
	MidCheckWidget* widget = new MidCheckWidget();
	auto dlg = GUIHelper::createDialog(widget, "MID clash detection");
	dlg->exec();
}

void MainWindow::on_actionVariantValidation_triggered()
{
	VariantValidationWidget* widget = new VariantValidationWidget();
	auto dlg = GUIHelper::createDialog(widget, "Variant validation");
	dlg->exec();
}

void MainWindow::on_actionChangePassword_triggered()
{
	PasswordDialog dlg(this);
	if(dlg.exec()==QDialog::Accepted)
	{
		NGSD db;
		db.setPassword(LoginManager::userId(), dlg.password());
	}
}

void MainWindow::on_actionStudy_triggered()
{
	DBTableAdministration* widget = new DBTableAdministration("study");
	auto dlg = GUIHelper::createDialog(widget, "Study administration");
	addModelessDialog(dlg);
}

void MainWindow::on_actionGaps_triggered()
{
	GapClosingDialog dlg(this);
	dlg.exec();
}

void MainWindow::on_actionReplicateNGSD_triggered()
{
	try
	{
		LoginManager::checkRoleIn(QStringList{"admin"});
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Permissions error", e.message());
		return;
	}

	NGSDReplicationWidget* widget = new NGSDReplicationWidget(this);
	auto dlg = GUIHelper::createDialog(widget, "Replicate NGSD (hg19 to hg38)");
	dlg->exec();
}

void MainWindow::on_actionPrepareGhgaUpload_triggered()
{
	GHGAUploadDialog dlg(this);
	dlg.exec();
}

void MainWindow::on_actionMaintenance_triggered()
{
	try
	{
		LoginManager::checkRoleIn(QStringList{"admin"});
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Permissions error", e.message());
		return;
	}


	MaintenanceDialog* dlg = new MaintenanceDialog(this);
	dlg->exec();
}

void MainWindow::on_actionNotifyUsers_triggered()
{
	try
	{
		LoginManager::checkRoleIn(QStringList{"admin"});
	}
	catch (Exception& e)
	{
		QMessageBox::warning(this, "Permissions error", e.message());
		return;
	}

	NGSD db;
	QStringList to, body;
	to = db.getValues("SELECT email from user WHERE user_role<>'special' AND active='1'");

	QStringList excluded_emails = {"restricted_user@med.uni-tuebingen.de"};
	foreach (QString email, excluded_emails)
	{
		int email_pos = to.indexOf(email);
		if (email_pos>-1) to.removeAt(email_pos);
	}

	QString subject = "GSvar update";
	body << "Dear all,";
	body << "";
	body << "";
	body << "Best regards,";
	body << LoginManager::userName();

	EmailDialog dlg(this, to, subject, body);
	dlg.exec();
}

void MainWindow::on_actionCohortAnalysis_triggered()
{
	CohortAnalysisWidget* widget = new CohortAnalysisWidget(this);
	auto dlg = GUIHelper::createDialog(widget, "Cohort analysis");
	addModelessDialog(dlg);
}



void MainWindow::on_actionGenderXY_triggered()
{
	ExternalToolDialog dialog("Determine gender", "xy", this);
	dialog.exec();
}

void MainWindow::on_actionGenderHet_triggered()
{
	ExternalToolDialog dialog("Determine gender", "hetx", this);
	dialog.exec();
}

void MainWindow::on_actionGenderSRY_triggered()
{
	ExternalToolDialog dialog("Determine gender", "sry", this);
	dialog.exec();
}

void MainWindow::on_actionStatisticsBED_triggered()
{
	ExternalToolDialog dialog("BED file information", "", this);
	dialog.exec();
}

void MainWindow::on_actionSampleSimilarityGSvar_triggered()
{
	ExternalToolDialog dialog("Sample similarity", "gsvar", this);
	dialog.exec();
}

void MainWindow::on_actionSampleSimilarityVCF_triggered()
{
	ExternalToolDialog dialog("Sample similarity", "vcf", this);
	dialog.exec();
}

void MainWindow::on_actionSampleSimilarityBAM_triggered()
{
	ExternalToolDialog dialog("Sample similarity", "bam", this);
	dialog.exec();
}

void MainWindow::on_actionSampleAncestry_triggered()
{
	ExternalToolDialog dialog("Sample ancestry", "", this);
	dialog.exec();
}

void MainWindow::on_actionAnalysisStatus_triggered()
{
	//check if already open
	for (int t=0; t<ui_.tabs->count(); ++t)
	{
		if (ui_.tabs->tabText(t)=="Analysis status")
		{
			ui_.tabs->setCurrentIndex(t);
			return;
		}
	}

	//open new
	AnalysisStatusWidget* widget = new AnalysisStatusWidget(this);
	connect(widget, SIGNAL(loadFile(QString)), this, SLOT(loadFile(QString)));
	openTab(QIcon(":/Icons/Server.png"), "Analysis status", widget);
}

void MainWindow::on_actionGapsLookup_triggered()
{
	if (filename_=="") return;

	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE) return;

	QString ps_name = selectProcessedSample();
	if (ps_name.isEmpty()) return;

	//check low-coverage file exists
	QStringList low_cov_files = GlobalServiceProvider::fileLocationProvider().getLowCoverageFiles(false).filterById(ps_name).asStringList();
	if (low_cov_files.isEmpty())
	{
		QMessageBox::warning(this, "Gap lookup", "No look-up of gaps is possible!\nCould not find a low-coverage file for sample " + ps_name + ".");
		return;
	}
	if (low_cov_files.count()>1) Log::warn( "Several gap files found for " + ps_name + ".");

	//get gene name from user
	QString gene = QInputDialog::getText(this, "Display gaps", "Gene:").trimmed();
	if (gene=="") return;

	//check if gene is in target region
	if (LoginManager::active())
	{
		NGSD db;
		QString ps_id = db.processedSampleId(germlineReportSample());
		if (ps_id!="")
		{
			int sys_id = db.getValue("SELECT processing_system_id FROM processed_sample WHERE id=:0", true, ps_id).toInt();
			BedFile sys_regions = GlobalServiceProvider::database().processingSystemRegions(sys_id, false);
			if (!sys_regions.isEmpty())
			{
				BedFile region = db.geneToRegions(gene.toUtf8(), Transcript::ENSEMBL, "gene");
				region.merge();
				if (region.count()==0)
				{
					QMessageBox::warning(this, "Precalculated gaps for gene", "Error:\nCould not convert gene symbol '" + gene + "' to a genomic region.\nIs this a HGNC-approved gene name with associated transcripts?");
					return;
				}

				region.intersect(sys_regions);
				if (region.count()==0)
				{
					QMessageBox::warning(this, "Precalculated gaps for gene", "Error:\nGene '" + gene + "' locus does not overlap with sample target region!");
					return;
				}
			}
		}
	}

	//look up data in report
	QStringList output;
	QStringList lines = Helper::loadTextFile(low_cov_files[0], true);
	foreach(QString line, lines)
	{
		QStringList parts = line.split('\t');
		if(parts.count()==4 && parts[3].contains(gene, Qt::CaseInsensitive))
		{
			double size_kb = (parts[2].toDouble() - parts[1].toDouble()) / 1000.0;
			output.append(line + "\t" + QString::number(size_kb, 'f', 3) + " kb");
		}
	}

	//show output
	QTextEdit* edit = new QTextEdit();
	edit->setText(output.join("\n"));
	edit->setMinimumWidth(500);
	edit->setWordWrapMode(QTextOption::NoWrap);
	edit->setReadOnly(true);
	auto dlg = GUIHelper::createDialog(edit, "Gaps of gene '" + gene + "' from low-coverage BED file for sample " + ps_name);
	dlg->exec();
}

void MainWindow::on_actionGapsRecalculate_triggered()
{
	if (filename_=="") return;

	//only available for gmerline and somatic single sample
	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE && type!=SOMATIC_SINGLESAMPLE) return;


	//check for BAM file
	QString ps = type==SOMATIC_SINGLESAMPLE ? variants_.getSampleHeader()[0].id : germlineReportSample();
	QStringList bams = GlobalServiceProvider::fileLocationProvider().getBamFiles(false).filterById(ps).asStringList();
	if (bams.empty())
	{
		QMessageBox::warning(this, "Gaps error", "No BAM file found for sample " + ps + "!");
		return;
	}

	//determine ROI name, ROI and gene list
	BedFile roi;
	GeneSet genes;

	//check for ROI file
	if (ui_.filters->targetRegion().isValid())
	{
		roi = ui_.filters->targetRegion().regions;
		genes = ui_.filters->targetRegion().genes;
	}
	else if (LoginManager::active())
	{
		QMessageBox::StandardButton btn = QMessageBox::information(this, "Gaps error", "No target region filter set!<br>Do you want to determine gaps for exon +- 20 bases of a specific gene?", QMessageBox::Yes, QMessageBox::No);
		if (btn!=QMessageBox::Yes) return;

		QByteArray symbol = selectGene().toUtf8();
		if (symbol=="") return;

		QApplication::setOverrideCursor(Qt::BusyCursor);

		roi = NGSD().geneToRegions(symbol, Transcript::ENSEMBL, "exon", true, false);
		roi.extend(20);
		roi.merge();

		genes << symbol;

		QApplication::restoreOverrideCursor();
	}
	else
	{
		QMessageBox::warning(this, "Gaps error", "No target region filter set!");
		return;
	}

	//show dialog
	QStringList low_covs = GlobalServiceProvider::fileLocationProvider().getLowCoverageFiles(false).filterById(ps).asStringList();
	low_covs << ""; //add empty string in case there is no low-coverage file > this case is handled inside the dialog
	if (low_covs.count()>1) Log::warn("Several gap files found for " + ps + ".");
	GapDialog dlg(this, ps, bams[0], low_covs[0], roi, genes);
	dlg.exec();
}

void MainWindow::exportVCF()
{
	try
	{
		QApplication::setOverrideCursor(Qt::BusyCursor);

		//generate GSvar with variants passing filter only
		VariantList selected_variants;
		selected_variants.copyMetaData(variants_);
		for(int i=0; i<variants_.count(); ++i)
		{
			if (!filter_result_.passing(i)) continue;
			selected_variants.append(variants_[i]);
		}

		//convert to VCF
		QString ref_genome = Settings::string("reference_genome", false);
		VcfFile vcf_file = VcfFile::fromGSvar(selected_variants, ref_genome);

		//store
		QString folder = Settings::path("gsvar_variant_export_folder", true);
		QString file_name = folder + QDir::separator() + QFileInfo(filename_).baseName() + "_export_" + QDate::currentDate().toString("yyyyMMdd") + "_" + Helper::userName() + ".vcf";

		file_name = QFileDialog::getSaveFileName(this, "Export VCF", file_name, "VCF (*.vcf);;All files (*.*)");
		if (file_name!="")
		{
			vcf_file.store(file_name);
			QApplication::restoreOverrideCursor();
			QMessageBox::information(this, "VCF export", "Exported VCF file with " + QString::number(vcf_file.count()) + " variants.");
		}
		else
		{
			QApplication::restoreOverrideCursor();
		}
	}
	catch(Exception& e)
	{
		QApplication::restoreOverrideCursor();

		QMessageBox::warning(this, "VCF export error", e.message());
	}
}

void MainWindow::exportGSvar()
{
	try
	{
		QApplication::setOverrideCursor(Qt::BusyCursor);

		//create new GSvar file with passing variants
		VariantList output;
		output.copyMetaData(variants_);
		for(int i=0; i<variants_.count(); ++i)
		{
			if (filter_result_.passing(i))
			{
				output.append(variants_[i]);
			}
		}

		//store
		QString folder = Settings::path("gsvar_variant_export_folder", true);
		QString file_name = folder + QDir::separator() + QFileInfo(filename_).baseName() + "_export_" + QDate::currentDate().toString("yyyyMMdd") + "_" + Helper::userName() + ".GSvar";

		file_name = QFileDialog::getSaveFileName(this, "Export GSvar", file_name, "VCF (*.GSvar);;All files (*.*)");
		if (file_name!="")
		{
			output.store(file_name);
			QApplication::restoreOverrideCursor();
			QMessageBox::information(this, "GSvar export", "Exported GSvar file with " + QString::number(output.count()) + " variants.");
		}
		else
		{
			QApplication::restoreOverrideCursor();
		}

	}
	catch(Exception& e)
	{
		QApplication::restoreOverrideCursor();

		QMessageBox::warning(this, "GSvar export error", e.message());
	}
}

void MainWindow::on_actionPreferredTranscripts_triggered()
{
	PreferredTranscriptsWidget* widget = new PreferredTranscriptsWidget();
	auto dlg = GUIHelper::createDialog(widget, "Preferred transcripts");
	dlg->exec();

	//re-load preferred transcripts from NGSD
	GSvarHelper::preferredTranscripts(true);
}

void MainWindow::on_actionEditSomaticGeneRoles_triggered()
{
	DBTableAdministration* table = new DBTableAdministration("somatic_gene_role");
	auto dlg = GUIHelper::createDialog(table, "Somatic Gene Roles");
	addModelessDialog(dlg);
}

void MainWindow::on_actionEditSomaticPathways_triggered()
{
	DBTableAdministration* table = new DBTableAdministration("somatic_pathway");
	auto dlg = GUIHelper::createDialog(table, "Somatic pathways");
	addModelessDialog(dlg);
}

void MainWindow::on_actionEditSomaticPathwayGeneAssociations_triggered()
{
	DBTableAdministration* table = new DBTableAdministration("somatic_pathway_gene");
	auto dlg = GUIHelper::createDialog(table, "Somatic pathways-gene associations");
	addModelessDialog(dlg);
}

void MainWindow::on_actionOpenDocumentation_triggered()
{
	QDesktopServices::openUrl(QUrl("https://github.com/imgag/ngs-bits/tree/master/doc/GSvar/index.md"));
}

void MainWindow::on_actionConvertHgnc_triggered()
{
	ApprovedGenesDialog dlg(this);
	dlg.exec();
}

void MainWindow::on_actionPhenoToGenes_triggered()
{
	try
	{
		PhenoToGenesDialog dlg(this);
		dlg.exec();
	}
	catch (DatabaseException& e)
	{
		QMessageBox::warning(this, "Database error", e.message());
	}
}

void MainWindow::on_actionGenesToRegions_triggered()
{
	GenesToRegionsDialog dlg(this);
	dlg.exec();
}

void MainWindow::openSubpanelDesignDialog(const GeneSet& genes)
{
	SubpanelDesignDialog dlg(this);
	dlg.setGenes(genes);

	dlg.exec();

	if (dlg.lastCreatedSubPanel()!="")
	{
		//update target region list
		ui_.filters->loadTargetRegions();

		//optionally use sub-panel as target regions
		if (QMessageBox::question(this, "Use sub-panel?", "Do you want to set the sub-panel as target region?")==QMessageBox::Yes)
		{
			ui_.filters->setTargetRegionByDisplayName(dlg.lastCreatedSubPanel());
		}
	}
}

void MainWindow::on_actionManageSubpanels_triggered()
{
	SubpanelArchiveDialog dlg(this);
	dlg.exec();
	if (dlg.changedSubpanels())
	{
		ui_.filters->loadTargetRegions();
	}
}

QString MainWindow::nobr()
{
	return "<p style='white-space:pre; margin:0; padding:0;'>";
}

void MainWindow::uploadToClinvar(int variant_index1, int variant_index2)
{
	if (!LoginManager::active()) return;

	try
	{
		if(variant_index1 < 0)
		{
			THROW(ArgumentException, "A valid variant index for the first variant has to be provided!");
		}
		//abort if API key is missing
		if(Settings::string("clinvar_api_key", true).trimmed().isEmpty())
		{
			THROW(ProgrammingException, "ClinVar API key is needed, but not found in settings.\nPlease inform the bioinformatics team");
		}

		NGSD db;

		//(1) prepare data as far as we can
		ClinvarUploadData data;
		data.processed_sample = germlineReportSample();
		data.variant_type1 = VariantType::SNVS_INDELS;
		if(variant_index2 < 0)
		{
			//Single variant submission
			data.submission_type = ClinvarSubmissionType::SingleVariant;
			data.variant_type2 = VariantType::INVALID;
		}
		else
		{
			//CompHet variant submission
			data.submission_type = ClinvarSubmissionType::CompoundHeterozygous;
			data.variant_type2 = VariantType::SNVS_INDELS;
		}

		QString sample_id = db.sampleId(data.processed_sample);
		SampleData sample_data = db.getSampleData(sample_id);


		//get disease info
		data.disease_info = db.getSampleDiseaseInfo(sample_id, "OMIM disease/phenotype identifier");
		data.disease_info.append(db.getSampleDiseaseInfo(sample_id, "Orpha number"));
		if (data.disease_info.length() < 1)
		{
			INFO(InformationMissingException, "The sample has to have at least one OMIM or Orphanet disease identifier to publish a variant in ClinVar.");
		}

		// get affected status
		data.affected_status = sample_data.disease_status;

		//get phenotype(s)
		data.phenos = sample_data.phenotypes;

		//get variant info
		data.snv1 = variants_[variant_index1];
		if(data.submission_type == ClinvarSubmissionType::CompoundHeterozygous) data.snv2 = variants_[variant_index2];

		// get report info
		if (!report_settings_.report_config.data()->exists(VariantType::SNVS_INDELS, variant_index1))
		{
			INFO(InformationMissingException, "The variant 1 has to be in the report configuration to be published!");
		}
		data.report_variant_config1 = report_settings_.report_config.data()->get(VariantType::SNVS_INDELS, variant_index1);
		if(data.submission_type == ClinvarSubmissionType::CompoundHeterozygous)
		{
			if (!report_settings_.report_config.data()->exists(VariantType::SNVS_INDELS, variant_index2))
			{
				INFO(InformationMissingException, "The variant 2 has to be in the report configuration to be published!");
			}
			data.report_variant_config2 = report_settings_.report_config.data()->get(VariantType::SNVS_INDELS, variant_index2);
		}



		//update classification
		data.report_variant_config1.classification = db.getClassification(data.snv1).classification;
		if (data.report_variant_config1.classification.trimmed().isEmpty() || (data.report_variant_config1.classification.trimmed() == "n/a"))
		{
			INFO(InformationMissingException, "The variant 1 has to be classified to be published!");
		}
		if(data.submission_type == ClinvarSubmissionType::CompoundHeterozygous)
		{
			data.report_variant_config2.classification = db.getClassification(data.snv2).classification;
			if (data.report_variant_config2.classification.trimmed().isEmpty() || (data.report_variant_config2.classification.trimmed() == "n/a"))
			{
				INFO(InformationMissingException, "The variant 2 has to be classified to be published!");
			}
		}

		//genes
		int gene_idx = variants_.annotationIndexByName("gene");
		data.genes = GeneSet::createFromText(data.snv1.annotations()[gene_idx], ',');
		if(data.submission_type == ClinvarSubmissionType::CompoundHeterozygous) data.genes <<  GeneSet::createFromText(data.snv2.annotations()[gene_idx], ',');

		//determine NGSD ids of variant and report variant for variant 1
		QString var_id = db.variantId(data.snv1, false);
		if (var_id == "")
		{
			INFO(InformationMissingException, "The variant 1 has to be in NGSD and part of a report config to be published!");
		}
		data.variant_id1 = Helper::toInt(var_id);
		//extract report variant id
		int rc_id = db.reportConfigId(db.processedSampleId(data.processed_sample));
		if (rc_id == -1 )
		{
			THROW(DatabaseException, "Could not determine report config id for sample " + data.processed_sample + "!");
		}

		data.report_variant_config_id1 = db.getValue("SELECT id FROM report_configuration_variant WHERE report_configuration_id=" + QString::number(rc_id) + " AND variant_id="
													 + QString::number(data.variant_id1), false).toInt();

		if(data.submission_type == ClinvarSubmissionType::CompoundHeterozygous)
		{
			//determine NGSD ids of variant and report variant for variant 2
			var_id = db.variantId(data.snv2, false);
			if (var_id == "")
			{
				INFO(InformationMissingException, "The variant 2 has to be in NGSD and part of a report config to be published!");
			}
			data.variant_id2 = Helper::toInt(var_id);

			//extract report variant id
			data.report_variant_config_id2 = db.getValue("SELECT id FROM report_configuration_variant WHERE report_configuration_id=" + QString::number(rc_id) + " AND variant_id="
														 + QString::number(data.variant_id2), false).toInt();
		}


		// (2) show dialog
		ClinvarUploadDialog dlg(this);
		dlg.setData(data);
		dlg.exec();
	}
	catch(Exception& e)
	{
		GUIHelper::showException(this, e, "ClinVar submission error");
	}
}

void MainWindow::dragEnterEvent(QDragEnterEvent* e)
{
	if (!e->mimeData()->hasFormat("text/uri-list")) return;
	if (e->mimeData()->urls().count()!=1) return;
	QUrl url = e->mimeData()->urls().at(0);
	if (!url.isLocalFile()) return;

	QString filename = url.toLocalFile();
	if (QFile::exists(filename) && filename.endsWith(".GSvar"))
	{
		e->acceptProposedAction();
	}
}

void MainWindow::dropEvent(QDropEvent* e)
{
	loadFile(e->mimeData()->urls().first().toLocalFile());
	e->accept();
}

void MainWindow::closeEvent(QCloseEvent* event)
{
	//unload the data
	loadFile();

    //close user session on the server
    if (ClientHelper::isClientServerMode()) performLogout();

	//here one could cancel closing the window by calling event->ignore()

	event->accept();
}

void MainWindow::refreshVariantTable(bool keep_widths)
{
	QApplication::setOverrideCursor(Qt::BusyCursor);

	QTime timer;
	timer.start();

	//apply filters
	applyFilters(false);
	int passing_variants = filter_result_.countPassing();
	QString status = QString::number(passing_variants) + " of " + QString::number(variants_.count()) + " variants passed filters.";
	int max_variants = 10000;
	if (passing_variants>max_variants)
	{
		status += " Displaying " + QString::number(max_variants) + " variants only!";
	}
	ui_.statusBar->showMessage(status);

	Log::perf("Applying all filters took ", timer);
	timer.start();

	//force update of variant details widget
	var_last_ = -1;

	//update variant table
	QList<int> col_widths = ui_.vars->columnWidths();
	AnalysisType type = variants_.type();
	if (type==SOMATIC_SINGLESAMPLE || type==SOMATIC_PAIR || type==CFDNA)
	{
		ui_.vars->update(variants_, filter_result_, somatic_report_settings_, max_variants);
	}
	else if (type==GERMLINE_SINGLESAMPLE || type==GERMLINE_TRIO || type==GERMLINE_MULTISAMPLE)
	{
		ui_.vars->update(variants_, filter_result_, report_settings_, max_variants);
	}
	else
	{
		THROW(ProgrammingException, "Unsupported analysis type in refreshVariantTable!");
	}

	ui_.vars->adaptRowHeights();
	if (keep_widths)
	{
		ui_.vars->setColumnWidths(col_widths);
	}
	else
	{
		ui_.vars->adaptColumnWidths();
	}
	QApplication::restoreOverrideCursor();

	Log::perf("Updating variant table took ", timer);
}


void MainWindow::varHeaderContextMenu(QPoint pos)
{
	if (!LoginManager::active()) return;

	//get variant index
	int index = ui_.vars->selectedVariantIndex();
	if(index==-1) return; //several variants selected

	//set up menu
	QMenu menu(ui_.vars->verticalHeader());
	QAction* a_edit = menu.addAction(QIcon(":/Icons/Report.png"), "Add/edit report configuration");
	QAction* a_delete =menu.addAction(QIcon(":/Icons/Remove.png"), "Delete report configuration");

	if(germlineReportSupported())
	{
		a_delete->setEnabled(!report_settings_.report_config->isFinalized() && report_settings_.report_config->exists(VariantType::SNVS_INDELS, index));
	}
	else if(somaticReportSupported())
	{
		 a_delete->setEnabled(somatic_report_settings_.report_config.exists(VariantType::SNVS_INDELS, index));
	}
	else
	{
		a_delete->setEnabled(false);
	}

	//exec menu
	pos = ui_.vars->verticalHeader()->viewport()->mapToGlobal(pos);
	QAction* action = menu.exec(pos);
	if (action==nullptr) return;

	//actions
	if (action==a_edit)
	{
		editVariantReportConfiguration(index);
	}
	else if (action==a_delete)
	{
		if(germlineReportSupported())
		{
			report_settings_.report_config->remove(VariantType::SNVS_INDELS, index);
		}
		else
		{
			somatic_report_settings_.report_config.remove(VariantType::SNVS_INDELS, index);
		}
		updateReportConfigHeaderIcon(index);
	}
}

void MainWindow::registerCustomContextMenuActions()
{
	bool  ngsd_user_logged_in = LoginManager::active();

	QList<QAction*> actions;
	context_menu_actions_.seperator = new QAction("---");

	//NGSD report configuration
	context_menu_actions_.a_report_edit = new QAction(QIcon(":/Icons/Report.png"), "Add/edit report configuration");
	context_menu_actions_.a_report_edit->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_report_edit;

	context_menu_actions_.a_report_del = new QAction(QIcon(":/Icons/Remove.png"), "Delete report configuration");
	context_menu_actions_.a_report_del->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_report_del;
	actions << context_menu_actions_.seperator;

	//NGSD variant options
	context_menu_actions_.a_var_class = new QAction("Edit classification");
	context_menu_actions_.a_var_class->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_var_class;

	context_menu_actions_.a_var_class_somatic = new QAction("Edit classification  (somatic)");
	context_menu_actions_.a_var_class_somatic->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_var_class_somatic;
	context_menu_actions_.a_var_interpretation_somatic = new QAction("Edit VICC interpretation (somatic)");
	context_menu_actions_.a_var_interpretation_somatic->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_var_interpretation_somatic;

	context_menu_actions_.a_var_comment = new QAction("Edit comment");
	context_menu_actions_.a_var_comment->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_var_comment;
	context_menu_actions_.a_var_val = new QAction("Perform variant validation");
	context_menu_actions_.a_var_val->setEnabled(ngsd_user_logged_in);
	actions << context_menu_actions_.a_var_val;
	actions << context_menu_actions_.seperator;

	ui_.vars->addCustomContextMenuActions(actions);
	connect(ui_.vars, SIGNAL(customActionTriggered(QAction*,int)), this, SLOT(execContextMenuAction(QAction*,int)));
}


void MainWindow::execContextMenuAction(QAction* action, int index)
{
	//perform actions
	if (action == context_menu_actions_.a_report_edit)
	{
		editVariantReportConfiguration(index);
	}
	else if (action == context_menu_actions_.a_report_del)
	{
		if ((!report_settings_.report_config->isFinalized() && report_settings_.report_config->exists(VariantType::SNVS_INDELS, index)) || somatic_report_settings_.report_config.exists(VariantType::SNVS_INDELS, index))
		{
			if(germlineReportSupported())
			{
				report_settings_.report_config->remove(VariantType::SNVS_INDELS, index);
			}
			else if(somaticReportSupported())
			{
				somatic_report_settings_.report_config.remove(VariantType::SNVS_INDELS, index);
				storeSomaticReportConfig();
			}

			updateReportConfigHeaderIcon(index);
		}
		else
		{
			QMessageBox::information(this, "Report configuration error", "This variant is not part of the report configuration. It can not be deleted from the report!");
		}
	}
	else if (action == context_menu_actions_.a_var_class)
	{
		editVariantClassification(variants_, index);
	}
	else if (action == context_menu_actions_.a_var_class_somatic)
	{
		editVariantClassification(variants_, index, true);
	}
	else if (action == context_menu_actions_.a_var_interpretation_somatic)
	{
		editSomaticVariantInterpretation(variants_, index);
	}
	else if (action == context_menu_actions_.a_var_comment)
	{
		editVariantComment(index);
	}
	else if (action == context_menu_actions_.a_var_val)
	{
		editVariantValidation(index);
	}
}

void MainWindow::openAlamut(QAction* action)
{
	//documentation of the alamut API:
	// - http://www.interactive-biosoftware.com/doc/alamut-visual/2.14/accessing.html
	// - http://www.interactive-biosoftware.com/doc/alamut-visual/2.11/Alamut-HTTP.html
	// - http://www.interactive-biosoftware.com/doc/alamut-visual/2.14/programmatic-access.html
	QStringList parts = action->text().split(" ");
	if (parts.count()>=1)
	{
		QString value = parts[0];
		if (value=="BAM")
		{
			QStringList bams = GlobalServiceProvider::fileLocationProvider().getBamFiles(false).filterById(germlineReportSample()).asStringList();
			if (bams.empty()) return;
			value = "BAM<" + bams[0];
		}

		try
		{
			QString host = Settings::string("alamut_host");
			QString institution = Settings::string("alamut_institution");
			QString apikey = Settings::string("alamut_apikey");
			HttpHandler(true).get(host+"/search?institution="+institution+"&apikey="+apikey+"&request="+value);
		}
		catch (Exception& e)
		{
			GUIHelper::showException(this, e, "Communication with Alamut failed!");
		}
	}
}

void MainWindow::showMatchingCnvsAndSvs(BedLine v_reg)
{
	try
	{
		//determine overlapping genes
		NGSD db;
		GeneSet genes = db.genesOverlapping(v_reg.chr(), v_reg.start(), v_reg.end());
		if (genes.isEmpty()) THROW(Exception, "Could not find a gene overlapping the variant region (gene transcripts are not padded by 5000 bases)!");

		//determine overlapping genes region
		BedFile regions = db.genesToRegions(genes, Transcript::ENSEMBL, "gene", true);
		regions.overlapping(v_reg); //sometimes genes have several loci due to duplicate gene names > exclude those
		regions.merge();

		//check target region
		if (regions.count()==0) THROW(Exception, "Could not determine a target region overlapping variant from genes: " + genes.join(", "));
		if (regions.count()>1) THROW(Exception, "Several target regions overlapping variant from genes: " + genes.join(", "));

		//create table
		TsvFile table;
		table.addHeader("type");
		table.addHeader("variant");
		table.addHeader("genotype");
		table.addHeader("details");

		//select CNVs
		{
			const Chromosome& chr = regions[0].chr();
			int start = regions[0].start();
			int end = regions[0].end();
			const QByteArrayList& headers = cnvs_.annotationHeaders();
			for (int i=0; i<cnvs_.count(); ++i)
			{
				const CopyNumberVariant& v = cnvs_[i];
				if (v.overlapsWith(chr, start, end))
				{
					int cn = v.copyNumber(headers);
					QStringList row;
					row << (cn<=1 ? "CNV - DEL" : "CNV - DUP");
					row << v.toString();
					row << "cn="+QString::number(cn);
					row << "size="+QString::number(v.size()/1000.0, 'f', 3) + "kb regions="+QString::number(v.regions());
					table.addRow(row);
				}
			}
		}

		//select SVs
		{
			QList<QByteArray> headers = svs_.annotationHeaders();
			for (int i=0; i<svs_.count(); ++i)
			{
				const BedpeLine& v = svs_[i];
				if (v.intersectsWith(regions))
				{
					QStringList row;
					row << ("SV - " + StructuralVariantTypeToString(v.type()));
					row << v.toString(false);
					row << "genotype="+v.genotypeHumanReadable(headers);
					row << "size="+QString::number(v.size()/1000.0, 'f', 3)+"kb";
					table.addRow(row);
				}
			}
		}

		//show table
		TsvTableWidget* widget = new TsvTableWidget(table, this);
		QSharedPointer<QDialog> dlg = GUIHelper::createDialog(widget, "CNVs and SVs matching " + v_reg.toString(true));
		dlg->exec();
	}
	catch(Exception& e)
	{
        GUIHelper::showException(this, e, "Showing matching CNVs and SVs failed!");
    }
}

void MainWindow::closeAndLogout()
{
    if (ClientHelper::isClientServerMode()) performLogout();
    close();
}

void MainWindow::displayIgvHistoryTable(QPoint /*pos*/)
{
	//check if already present > bring to front
	QList<IgvLogWidget*> dialogs = findChildren<IgvLogWidget*>();
	if (!dialogs.isEmpty())
	{
		QDialog* dlg = qobject_cast<QDialog*>(dialogs.at(0)->parent());
		dlg->raise();
		return;
	}

	//open new dialog
    IgvLogWidget* widget = new IgvLogWidget(this);
	auto dlg = GUIHelper::createDialog(widget, "IGV command history");
	addModelessDialog(dlg);
}

void MainWindow::changeIgvIconToActive()
{
    igv_history_label_->setPixmap(QPixmap(":/Icons/IGV_active.png"));
    igv_history_label_->setToolTip("IGV is currently processing a command");
}

void MainWindow::changeIgvIconToNormal()
{
	igv_history_label_->setPixmap(QPixmap(":/Icons/IGV.png"));
	igv_history_label_->setToolTip("IGV is idle at the moment");
}

void MainWindow::on_actionVirusDetection_triggered()
{
	//get virus file
	QString ps_tumor = variants_.mainSampleName();
	NGSD db;
	QString ps_tumor_id = db.processedSampleId(ps_tumor, false);
	FileLocation virus_file = GlobalServiceProvider::database().processedSamplePath(ps_tumor_id, PathType::VIRAL);

	//show widget
	VirusDetectionWidget* widget = new VirusDetectionWidget(virus_file.filename);
	auto dlg = GUIHelper::createDialog(widget, "Virus detection");
	addModelessDialog(dlg);
}

void MainWindow::on_actionBurdenTest_triggered()
{
	BurdenTestWidget* widget = new BurdenTestWidget(this);

	auto dlg = GUIHelper::createDialog(widget, "Gene-based burden test");
	addModelessDialog(dlg);
}


void MainWindow::editVariantClassification(VariantList& variants, int index, bool is_somatic)
{
	try
	{
		Variant& variant = variants[index];

		//execute dialog
		ClassificationDialog dlg(this, variant, is_somatic);
		if (dlg.exec()!=QDialog::Accepted) return;

		//update NGSD
		NGSD db;

		ClassificationInfo class_info = dlg.classificationInfo();
		if(is_somatic)
		{
			db.setSomaticClassification(variant, class_info);

			//update variant list classification
			int i_som_class = variants.annotationIndexByName("somatic_classification");
			QString new_class = class_info.classification.replace("n/a", "");
			variant.annotations()[i_som_class] = new_class.toUtf8();

			markVariantListChanged(variant, "somatic_classification", new_class);

			//update variant list classification comment
			int i_som_class_comment = variants.annotationIndexByName("somatic_classification_comment");
			variant.annotations()[i_som_class_comment] = class_info.comments.toUtf8();

			markVariantListChanged(variant, "somatic_classification_comment", class_info.comments);

		}
		else //germline variants
		{
			db.setClassification(variant, variants_, class_info);

			//update variant list classification
			int i_class = variants.annotationIndexByName("classification");
			QString new_class = class_info.classification.replace("n/a", "");
			variant.annotations()[i_class] = new_class.toUtf8();

			markVariantListChanged(variant, "classification", new_class);

			//update variant list classification comment
			int i_class_comment = variants.annotationIndexByName("classification_comment");
			variant.annotations()[i_class_comment] = class_info.comments.toUtf8();

			markVariantListChanged(variant, "classification_comment", class_info.comments);

			//check if already uploaded to ClinVar
			QString var_id = db.variantId(variant);
			QString sample_id = db.sampleId(germlineReportSample());
			QString clinvar_class = db.getValue("SELECT class FROM  variant_publication WHERE variant_table='variant' AND db='ClinVar' AND sample_id='" + sample_id + "' AND variant_id='" + var_id + "' ORDER BY id DESC LIMIT 1").toString();
			if(!clinvar_class.isEmpty() && clinvar_class!=new_class)
			{
				//update on ClinVar
				int return_value = QMessageBox::information(this, "Clinvar upload required!", "Variant already uploaded to ClinVar. You should also update the classification there!", QMessageBox::Ok, QMessageBox::NoButton);
				if(return_value == QMessageBox::Ok)	uploadToClinvar(index);
			}
		}

		//update details widget and filtering
		ui_.variant_details->updateVariant(variants, index);
		refreshVariantTable();

	}
	catch (DatabaseException& e)
	{
		GUIHelper::showMessage("NGSD error", e.message());
		return;
	}
}

void MainWindow::editSomaticVariantInterpretation(const VariantList &vl, int index)
{
	SomaticVariantInterpreterWidget* interpreter = new SomaticVariantInterpreterWidget(index, vl, this);
	auto dlg = GUIHelper::createDialog(interpreter, "Somatic Variant Interpretation");
	connect(interpreter, SIGNAL(stored(int, QString, QString)), this, SLOT(updateSomaticVariantInterpretationAnno(int, QString, QString)) );

	dlg->exec();
}

void MainWindow::updateSomaticVariantInterpretationAnno(int index, QString vicc_interpretation, QString vicc_comment)
{
	int i_vicc = variants_.annotationIndexByName("NGSD_som_vicc_interpretation");
	variants_[index].annotations()[i_vicc] = vicc_interpretation.toUtf8();

	markVariantListChanged(variants_[index], "NGSD_som_vicc_interpretation", vicc_interpretation);

	int i_vicc_comment = variants_.annotationIndexByName("NGSD_som_vicc_comment");
	variants_[index].annotations()[i_vicc_comment] = vicc_comment.toUtf8();

	markVariantListChanged(variants_[index], "NGSD_som_vicc_comment", vicc_comment);

	//update details widget and filtering
	ui_.variant_details->updateVariant(variants_, index);
	refreshVariantTable();
}

void MainWindow::on_actionAnnotateSomaticVariantInterpretation_triggered()
{
	if (filename_.isEmpty()) return;
	if (!LoginManager::active()) return;
	AnalysisType type = variants_.type();
	if (type!=SOMATIC_SINGLESAMPLE && type!=SOMATIC_PAIR) return;

	int i_vicc = variants_.annotationIndexByName("NGSD_som_vicc_interpretation");
	int i_vicc_comment = variants_.annotationIndexByName("NGSD_som_vicc_comment");

	NGSD db;
	for(int i=0; i<variants_.count(); ++i)
	{
		//skip variants without VICC infos in NGSD
		SomaticViccData vicc_data = db.getSomaticViccData(variants_[i], false);
		if (vicc_data.created_by.isEmpty()) continue;

		//update score
		QByteArray vicc_score = SomaticVariantInterpreter::viccScoreAsString(vicc_data).toUtf8();
		if (vicc_score!=variants_[i].annotations()[i_vicc])
		{
			variants_[i].annotations()[i_vicc] = vicc_score;
			markVariantListChanged(variants_[i], "NGSD_som_vicc_interpretation", vicc_score);
		}

		//update comment
		QByteArray vicc_comment = vicc_data.comment.toUtf8();
		if (variants_[i].annotations()[i_vicc_comment]!=vicc_comment)
		{
			variants_[i].annotations()[i_vicc_comment]= vicc_comment;
			markVariantListChanged(variants_[i], "NGSD_som_vicc_comment", vicc_comment);
		}
	}

	//update details widget and filtering
	refreshVariantTable();
}

bool MainWindow::germlineReportSupported(bool require_ngsd)
{
	//no file loaded
	if (filename_.isEmpty()) return false;

	//user has to be logged in
	if (require_ngsd && !LoginManager::active()) return false;

	//single, trio or multi only
	AnalysisType type = variants_.type();
	if (type!=GERMLINE_SINGLESAMPLE && type!=GERMLINE_TRIO && type!=GERMLINE_MULTISAMPLE) return false;

	//multi-sample only with at least one affected
	if (type==GERMLINE_MULTISAMPLE && variants_.getSampleHeader().sampleColumns(true).count()<1) return false;

	//affected samples are in NGSD
	if (require_ngsd)
	{
		NGSD db;
		foreach(const SampleInfo& info, variants_.getSampleHeader())
		{
			if(info.isAffected())
			{
				if (db.processedSampleId(info.id.trimmed(), false)=="") return false;
			}
		}
	}

	return true;
}

QString MainWindow::germlineReportSample()
{
	if (!germlineReportSupported(false))
	{
		THROW(ProgrammingException, "germlineReportSample() cannot be used if germline report is not supported!");
	}

	//set sample for report
	while (germline_report_ps_.isEmpty())
	{
		//determine affected sample names
		QStringList affected_ps;
		foreach(const SampleInfo& info, variants_.getSampleHeader())
		{
			if(info.isAffected())
			{
				affected_ps << info.id.trimmed();
			}
		}

		if (affected_ps.isEmpty()) //no affected => error
		{
			THROW(ProgrammingException, "germlineReportSample() cannot be used if there is no affected sample!");
		}
		else if (affected_ps.count()==1) //one affected => auto-select
		{
			germline_report_ps_ = affected_ps[0];
		}
		else //several affected => let user select
		{
			bool ok = false;
			QString selected = QInputDialog::getItem(this, "Report sample", "processed sample used for report:", affected_ps, 0, false, &ok);
			if (ok)
			{
				germline_report_ps_ = selected;
			}
		}
	}

	return germline_report_ps_;
}

bool MainWindow::somaticReportSupported()
{
	return variants_.type()==SOMATIC_PAIR;
}

bool MainWindow::tumoronlyReportSupported()
{
	return variants_.type()==SOMATIC_SINGLESAMPLE;
}

void MainWindow::updateVariantDetails()
{
	int var_current = ui_.vars->selectedVariantIndex();
	if (var_current==-1) //no several variant => clear
	{
		ui_.variant_details->clear();
	}
	else if (var_current!=var_last_) //update variant details (if changed)
	{
		ui_.variant_details->updateVariant(variants_, var_current);
	}

	var_last_ = var_current;
}

void MainWindow::editVariantReportConfiguration(int index)
{
	if (!germlineReportSupported() && !somaticReportSupported())
	{
		QMessageBox::information(this, "Report configuration error", "Report configuration not supported for this type of analysis!");
		return;
	}

	NGSD db;

	if(germlineReportSupported()) //germline report configuration
	{
		//init/get config
		ReportVariantConfiguration var_config;
		bool report_settings_exist = report_settings_.report_config->exists(VariantType::SNVS_INDELS, index);
		if (report_settings_exist)
		{
			var_config = report_settings_.report_config->get(VariantType::SNVS_INDELS, index);
		}
		else
		{
			var_config.variant_index = index;
		}

		//get inheritance mode by gene
		const Variant& variant = variants_[index];
		QList<KeyValuePair> inheritance_by_gene;
		int i_genes = variants_.annotationIndexByName("gene", true, false);

		if (i_genes!=-1)
		{
			GeneSet genes = GeneSet::createFromText(variant.annotations()[i_genes], ',');
			foreach(const QByteArray& gene, genes)
			{
				GeneInfo gene_info = db.geneInfo(gene);
				inheritance_by_gene << KeyValuePair{gene, gene_info.inheritance};
			}
		}

		//exec dialog
		ReportVariantDialog dlg(variant.toString(false), inheritance_by_gene, var_config, this);
		dlg.setEnabled(!report_settings_.report_config->isFinalized());
		if (dlg.exec()!=QDialog::Accepted) return;


		//update config, GUI and NGSD
		report_settings_.report_config->set(var_config);
		updateReportConfigHeaderIcon(index);

		//force classification of causal variants
		if(var_config.causal)
		{
			const Variant& variant = variants_[index];
			ClassificationInfo classification_info = db.getClassification(variant);
			if (classification_info.classification=="" || classification_info.classification=="n/a")
			{
				QMessageBox::warning(this, "Variant classification required!", "Causal variants should have a classification!", QMessageBox::Ok, QMessageBox::NoButton);
				editVariantClassification(variants_, index);
			}

			//enforce ClinVar upload of class 4/5 variants
			classification_info = db.getClassification(variant);
			if (classification_info.classification=="4" || classification_info.classification=="5")
			{
				QList<int> publication_ids = db.getValuesInt("SELECT id FROM variant_publication WHERE variant_id='" + db.variantId(variants_[index]) + "'");
				if (publication_ids.isEmpty())
				{
					QMessageBox::information(this, "Clinvar upload required!", "Class 4 or 5 variants should be uploaded to ClinVar!", QMessageBox::Ok, QMessageBox::NoButton);
					uploadToClinvar(index);
				}
			}
		}
	}
	else if(somaticReportSupported()) //somatic report variant configuration
	{
		SomaticReportVariantConfiguration var_config;
		bool settings_exists = somatic_report_settings_.report_config.exists(VariantType::SNVS_INDELS, index);
		if(settings_exists)
		{
			var_config = somatic_report_settings_.report_config.get(VariantType::SNVS_INDELS, index);
		}
		else
		{
			var_config.variant_index = index;
		}

		SomaticReportVariantDialog* dlg = new SomaticReportVariantDialog(variants_[index].toString(), var_config, this);

		if(dlg->exec() != QDialog::Accepted) return;
		somatic_report_settings_.report_config.set(var_config);

		storeSomaticReportConfig();
		updateReportConfigHeaderIcon(index);
	}
}

void MainWindow::updateReportConfigHeaderIcon(int index)
{
	if(germlineReportSupported())
	{
		//report config-based filter is on => update whole variant list
		if (ui_.filters->reportConfigurationFilter()!=ReportConfigFilter::NONE)
		{
			refreshVariantTable();
		}
		else //no filter => refresh icon only
		{
			ui_.vars->updateVariantHeaderIcon(report_settings_, index);
		}
	}
	else if(somaticReportSupported())
	{
		if(!ui_.filters->targetRegion().isValid() || ui_.filters->filters().count() > 0)
		{
			refreshVariantTable();
		}
		else
		{
			ui_.vars->updateVariantHeaderIcon(somatic_report_settings_, index);
		}
	}
}

void MainWindow::markVariantListChanged(const Variant& variant, QString column, QString text)
{
	variants_changed_ << VariantListChange{variant, column, text};
}

void MainWindow::storeCurrentVariantList()
{
	QApplication::setOverrideCursor(Qt::BusyCursor);

	if (GlobalServiceProvider::fileLocationProvider().isLocal())
	{
		try
		{
			//store to temporary file
			QString tmp = filename_ + ".tmp";
			variants_.store(tmp);

			//copy temp
			QFile::remove(filename_);
			QFile::rename(tmp, filename_);

			variants_changed_.clear();
		}
		catch(Exception& e)
		{
			QApplication::restoreOverrideCursor();
			QMessageBox::warning(this, "Error storing GSvar file", "The GSvar file could not be stored:\n" + e.message());
		}
	}
	else
	{
		QJsonDocument json_doc = QJsonDocument();
		QJsonArray json_array;
		QJsonObject json_object;

		foreach(const VariantListChange& variant_changed, variants_changed_)
		{
			try
			{
				json_object.insert("variant", variant_changed.variant.toString());
				json_object.insert("column", variant_changed.column);
				json_object.insert("text", variant_changed.text);
				json_array.append(json_object);
			}
			catch (Exception& e)
			{
				QMessageBox::warning(this, "Could not process the changes to be sent to the server:", e.message());
			}
		}

		json_doc.setArray(json_array);

		QString ps_url_id;
		QList<QString> filename_parts = filename_.split("/");
		if (filename_parts.size()>3)
		{
			ps_url_id = filename_parts[filename_parts.size()-2];
		}

		try
		{
			HttpHeaders add_headers;
			add_headers.insert("Accept", "application/json");
			add_headers.insert("Content-Type", "application/json");
			add_headers.insert("Content-Length", QByteArray::number(json_doc.toJson().count()));

			QString reply = HttpHandler(true).put(
						ClientHelper::serverApiUrl() + "project_file?ps_url_id=" + ps_url_id + "&token=" + LoginManager::userToken(),
						json_doc.toJson(),
						add_headers
					);
		}
		catch (Exception& e)
		{
			QMessageBox::warning(this, "Could not reach the server:", e.message());
		}
	}

	QApplication::restoreOverrideCursor();
}

void MainWindow::checkPendingVariantValidations()
{
	if (!LoginManager::active()) return;

	NGSD db;
	QStringList vv_pending = db.getValues("SELECT id FROM variant_validation WHERE status='for reporting' AND user_id='" + LoginManager::userIdAsString() + "'");
	if (vv_pending.isEmpty()) return;

	showNotification("Variant validation: " + QString::number(vv_pending.count()) + " pending variants 'for reporting'!");
}

void MainWindow::showNotification(QString text)
{
	text = text.trimmed();

	//update tooltip
	QStringList tooltips = notification_label_->toolTip().split("\n", QString::SkipEmptyParts);
	if (!tooltips.contains(text)) tooltips.prepend(text);
	notification_label_->setToolTip(tooltips.join("<br>"));

	//show popup
	notification_label_->show();
	QPoint pos = ui_.statusBar->mapToGlobal(notification_label_->pos()) + QPoint(8,8);
	QToolTip::showText(pos, text);
}

void MainWindow::variantRanking()
{
	if (filename_.isEmpty()) return;
	if (!LoginManager::active()) return;

	QApplication::setOverrideCursor(Qt::BusyCursor);
	QString ps_name = germlineReportSample();
	try
	{
		NGSD db;

		//create phenotype list
		QHash<Phenotype, BedFile> phenotype_rois;
		QString sample_id = db.sampleId(ps_name);
		PhenotypeList phenotypes = ui_.filters->phenotypes();
		if (phenotypes.isEmpty())
		{
			phenotypes = db.getSampleData(sample_id).phenotypes;
		}
		foreach(const Phenotype& pheno, phenotypes)
		{
			//pheno > genes
			GeneSet genes = db.phenotypeToGenes(db.phenotypeIdByAccession(pheno.accession()), true);

			//genes > roi
			BedFile roi;
			foreach(const QByteArray& gene, genes)
			{
				if (!gene2region_cache_.contains(gene))
				{
					BedFile tmp = db.geneToRegions(gene, Transcript::ENSEMBL, "gene", true);
					tmp.clearAnnotations();
					tmp.merge();
					gene2region_cache_[gene] = tmp;
				}
				roi.add(gene2region_cache_[gene]);
			}
			roi.merge();

			phenotype_rois[pheno] = roi;
		}

		//score
		VariantScores::Parameters parameters;
		QString algorithm = sender()->objectName();
		VariantScores::Result result = VariantScores::score(algorithm, variants_, phenotype_rois, parameters);

		//update variant list
		VariantScores::annotate(variants_, result, true);
		ui_.filters->reset(true);
		ui_.filters->setFilter("GSvar score/rank");

		QApplication::restoreOverrideCursor();

		//show warnings
		if (result.warnings.count()>0)
		{
			QMessageBox::warning(this, "Variant ranking", "Please note the following warnings:\n" + result.warnings.join("\n"));
		}
	}
	catch(Exception& e)
	{
		QApplication::restoreOverrideCursor();
		QMessageBox::warning(this, "Ranking variants", "An error occurred:\n" + e.message());
	}
}

void MainWindow::clearSomaticReportSettings(QString ps_id_in_other_widget)
{
	if(!LoginManager::active()) return;
	QString this_ps_id = NGSD().processedSampleId(variants_.mainSampleName(),false);

	if(this_ps_id == "") return;

	if(this_ps_id != ps_id_in_other_widget) return; //skip if ps id of file is different than in other widget
	somatic_report_settings_ = SomaticReportSettings();
	refreshVariantTable();
}

void MainWindow::applyFilters(bool debug_time)
{
	try
	{
		//apply main filter
		QTime timer;
		timer.start();

		const FilterCascade& filter_cascade = ui_.filters->filters();

		filter_result_ = filter_cascade.apply(variants_, false, debug_time);

		ui_.filters->markFailedFilters();

		if (debug_time)
		{
			Log::perf("Applying annotation filters took ", timer);
			timer.start();
		}

		//roi filter
		if (ui_.filters->targetRegion().isValid())
		{
			FilterRegions::apply(variants_, ui_.filters->targetRegion().regions, filter_result_);

			if (debug_time)
			{
				Log::perf("Applying target region filter took ", timer);
				timer.start();
			}
		}

		//gene filter
		GeneSet genes_filter = ui_.filters->genes();
		if (!genes_filter.isEmpty())
		{
			FilterGenes filter;
			filter.setStringList("genes", genes_filter.toStringList());
			filter.apply(variants_, filter_result_);

			if (debug_time)
			{
				Log::perf("Applying gene filter took ", timer);
				timer.start();
			}
		}

		//text filter
		QByteArray text = ui_.filters->text();
		if (!text.isEmpty())
		{
			FilterAnnotationText filter;
			filter.setString("term", text);
			filter.setString("action", "FILTER");
			filter.apply(variants_, filter_result_);

			if (debug_time)
			{

				Log::perf("Applying text filter took ", timer);
				timer.start();
			}
		}

		//region filter
		QString region_text = ui_.filters->region();
		BedLine region = BedLine::fromString(region_text);
		if (!region.isValid()) //check if valid chr
		{
			Chromosome chr(region_text);
			if (chr.isNonSpecial())
			{
				region.setChr(chr);
				region.setStart(1);
				region.setEnd(999999999);
			}
		}
		if (region.isValid()) //valid region (chr,start, end or only chr)
		{
			BedFile tmp;
			tmp.append(region);
			FilterRegions::apply(variants_, tmp, filter_result_);

			if (debug_time)
			{
				Log::perf("Applying region filter took ", timer);
				timer.start();
			}
		}
		//phenotype selection changed => update ROI
		const PhenotypeList& phenos = ui_.filters->phenotypes();

		//update phenotypes for variant context menu search
		ui_.vars->updateActivePhenotypes(phenos);

		const PhenotypeSettings& pheno_settings = ui_.filters->phenotypeSettings();
		if (phenos!=last_phenos_ || pheno_settings!=last_pheno_settings_)
		{
			last_phenos_ = phenos;
			last_pheno_settings_ = pheno_settings;

			//convert phenotypes to genes
			NGSD db;
			GeneSet pheno_genes;
			int i = 0;
			foreach(const Phenotype& pheno, phenos)
			{
				GeneSet genes = db.phenotypeToGenesbySourceAndEvidence(db.phenotypeIdByAccession(pheno.accession()), pheno_settings.sources, pheno_settings.evidence_levels, true, false);

				if (pheno_settings.mode==PhenotypeCombimnationMode::MERGE || (pheno_settings.mode==PhenotypeCombimnationMode::INTERSECT && i==0))
				{
					pheno_genes << genes;
				}
				else
				{
					pheno_genes = pheno_genes.intersect(genes);
				}
				++i;
			}

			//convert genes to ROI (using a cache to speed up repeating queries)
			phenotype_roi_.clear();
			foreach(const QByteArray& gene, pheno_genes)
			{
				if (!gene2region_cache_.contains(gene))
				{
					BedFile tmp = db.geneToRegions(gene, Transcript::ENSEMBL, "gene", true);
					tmp.clearAnnotations();
					tmp.extend(5000);
					tmp.merge();
					gene2region_cache_[gene] = tmp;
				}
				phenotype_roi_.add(gene2region_cache_[gene]);
			}
			phenotype_roi_.merge();

			if (debug_time)
			{
				Log::perf("Updating phenotype filter took ", timer);
				timer.start();
			}
		}

		//phenotype filter
		if (!last_phenos_.isEmpty())
		{
			FilterRegions::apply(variants_, phenotype_roi_, filter_result_);

			if (debug_time)
			{
				Log::perf("Applying phenotype filter took ", timer);
				timer.start();
			}
		}

		//report configuration filter (show only variants with report configuration)
		ReportConfigFilter rc_filter = ui_.filters->reportConfigurationFilter();
		if (germlineReportSupported() && rc_filter!=ReportConfigFilter::NONE)
		{
			QSet<int> report_variant_indices = report_settings_.report_config->variantIndices(VariantType::SNVS_INDELS, false).toSet();
			for(int i=0; i<variants_.count(); ++i)
			{
				if (!filter_result_.flags()[i]) continue;

				if (rc_filter==ReportConfigFilter::HAS_RC)
				{
					filter_result_.flags()[i] = report_variant_indices.contains(i);
				}
				else if (rc_filter==ReportConfigFilter::NO_RC)
				{
					filter_result_.flags()[i] = !report_variant_indices.contains(i);
				}
			}
		}
		else if( somaticReportSupported() && rc_filter != ReportConfigFilter::NONE) //somatic report configuration filter (show only variants with report configuration)
		{
			QSet<int> report_variant_indices = somatic_report_settings_.report_config.variantIndices(VariantType::SNVS_INDELS, false).toSet();
			for(int i=0; i<variants_.count(); ++i)
			{
				if ( !filter_result_.flags()[i] ) continue;

				if (rc_filter==ReportConfigFilter::HAS_RC)
				{
					filter_result_.flags()[i] = report_variant_indices.contains(i);
				}
				else if (rc_filter==ReportConfigFilter::NO_RC)
				{
					filter_result_.flags()[i] = !report_variant_indices.contains(i);
				}
			}
		}

		//keep somatic variants that are marked with "include" in report settings (overrides possible filtering for that variant)
		if( somaticReportSupported() && rc_filter != ReportConfigFilter::NO_RC)
		{
			foreach(int index, somatic_report_settings_.report_config.variantIndices(VariantType::SNVS_INDELS, false))
			{
				filter_result_.flags()[index] = filter_result_.flags()[index] || somatic_report_settings_.report_config.variantConfig(index, VariantType::SNVS_INDELS).showInReport();
			}
		}
	}
	catch(Exception& e)
	{
		QMessageBox::warning(this, "Filtering error", e.message() + "\nA possible reason for this error is an outdated variant list.\nPlease re-run the annotation steps for the analysis!");

		filter_result_ = FilterResult(variants_.count(), false);
	}
}

void MainWindow::addToRecentSamples(QString ps)
{
	//update settings
	QStringList recent_samples = Settings::stringList("recent_samples", true);
	recent_samples.removeAll(ps);
	recent_samples.prepend(ps);
	recent_samples = recent_samples.mid(0,10);

	Settings::setStringList("recent_samples", recent_samples);

	//update GUI
	updateRecentSampleMenu();
}


void MainWindow::updateRecentSampleMenu()
{
	QStringList recent_samples = Settings::stringList("recent_samples", true);

	QMenu* menu = new QMenu();
	foreach(const QString& sample, recent_samples)
	{
		menu->addAction(sample, this, SLOT(openRecentSample()));
	}
	ui_.actionRecent->setMenu(menu);
}

void MainWindow::updateIGVMenu()
{
	QStringList entries = Settings::stringList("igv_menu");
	if (entries.count()==0)
	{
		ui_.menuOpenCustomTrack->addAction("No custom entries in INI file!");
	}
	else
	{
		foreach(QString entry, entries)
		{
			QStringList parts = entry.trimmed().split("\t");
			if(parts.count()!=3) continue;

			QString filename = parts[2];

			//add to menu "open custom track"
			QAction* action = ui_.menuOpenCustomTrack->addAction(parts[0], this, SLOT(openCustomIgvTrack()));
			action->setToolTip(filename); //file path is taken from tooltip when opening track
			if (!QFile::exists(filename))
			{
				action->setEnabled(false);
				action->setText(action->text() + " (missing)");
			}
		}
	}
}

void MainWindow::updateNGSDSupport()
{
	//init
	bool ngsd_user_logged_in = LoginManager::active();

	//toolbar
	ui_.actionOpenProcessedSampleTabByName->setEnabled(ngsd_user_logged_in);
	ui_.actionOpenSequencingRunTabByName->setEnabled(ngsd_user_logged_in);
	ui_.actionOpenGeneTabByName->setEnabled(ngsd_user_logged_in);
	ui_.actionOpenVariantTab->setEnabled(ngsd_user_logged_in);
	ui_.actionOpenProjectTab->setEnabled(ngsd_user_logged_in);
	ui_.actionOpenProcessingSystemTab->setEnabled(ngsd_user_logged_in);
	ui_.report_btn->setEnabled(ngsd_user_logged_in);
	ui_.actionAnalysisStatus->setEnabled(ngsd_user_logged_in);
	ui_.actionReanalyze->setEnabled(ngsd_user_logged_in);
	ui_.actionGapsRecalculate->setEnabled(ngsd_user_logged_in);
	ui_.actionGeneSelector->setEnabled(ngsd_user_logged_in);
	ui_.actionSampleSearch->setEnabled(ngsd_user_logged_in);
	ui_.actionRunOverview->setEnabled(ngsd_user_logged_in);
	ui_.actionConvertHgvsToGSvar->setEnabled(ngsd_user_logged_in);
	ui_.actionRegionToGenes->setEnabled(ngsd_user_logged_in);
	ui_.actionGapsRecalculate->setEnabled(ngsd_user_logged_in);
	ui_.actionAnnotateSomaticVariantInterpretation->setEnabled(ngsd_user_logged_in);

	//NGSD menu
	ui_.menuNGSD->setEnabled(ngsd_user_logged_in);
	ui_.actionDesignSubpanel->setEnabled(ngsd_user_logged_in);

	//other actions
	ui_.actionOpenByName->setEnabled(ngsd_user_logged_in);
	ui_.ps_details->setEnabled(ngsd_user_logged_in);
	ui_.vars_ranking->setEnabled(ngsd_user_logged_in);

	ui_.filters->updateNGSDSupport();

	//disable certain actions/buttons for restricted users
	if (ngsd_user_logged_in)
	{
		NGSD db;
		if (db.userRoleIn(LoginManager::userLogin(), QStringList{"user_restricted"}))
		{
			auto actions = ui_.menuAdmin->actions();
			foreach(QAction* action, actions)
			{
				if (action!=ui_.actionChangePassword)
				{
					action->setEnabled(false);
				}
			}
		}
	}
}

void MainWindow::openRecentSample()
{
	QAction* action = qobject_cast<QAction*>(sender());
	openProcessedSampleFromNGSD(action->text());
}

QString MainWindow::normalSampleName()
{
	if(variants_.type() != AnalysisType::SOMATIC_PAIR) return "";

	foreach(const SampleInfo& info, variants_.getSampleHeader())
	{
		if (!info.isTumor()) return info.id;
	}

	return "";
}

QString MainWindow::getFileSelectionItem(QString window_title, QString label_text, QStringList file_list, bool *ok)
{
	QStringList rna_count_files_displayed_names = file_list;
	if (ClientHelper::isClientServerMode())
	{
		rna_count_files_displayed_names.clear();
		foreach (QString full_name, file_list)
		{
			rna_count_files_displayed_names << QUrl(full_name).fileName();
		}
	}

	QString selected_file_name = QInputDialog::getItem(this, window_title, label_text, rna_count_files_displayed_names, 0, false, ok);
	if (ClientHelper::isClientServerMode())
	{
		int selection_index = rna_count_files_displayed_names.indexOf(selected_file_name);
		if (selection_index>-1) return file_list[selection_index];
		QMessageBox::warning(this, "File URL not found", "Could not find the URL for the selected file!");
		return "";
	}

    return selected_file_name;
}

void MainWindow::performLogout()
{
    if (!LoginManager::active()) return;

    LoginManager::logout();
}



