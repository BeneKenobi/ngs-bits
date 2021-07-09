#include <QCoreApplication>
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QCommandLineParser>
#include "HttpsServer.h"
#include "ServerHelper.h"
#include "EndpointController.h"
#include "EndpointHandler.h"

int log_level = 3;

QFile gsvar_server_log_file("gsvar-server-log.txt");

void interceptLogMessage(QtMsgType type, const QMessageLogContext &, const QString &msg)
{
	QString time_stamp = QDate::currentDate().toString("dd/MM/yyyy") + " " + QTime::currentTime().toString("hh:mm:ss:zzz");
	QString log_statement = "";
	int msg_level = 0;
	switch (type) {
		case QtCriticalMsg:
			msg_level = 0;
			log_statement = QString("%1 - [Critical] %2").arg(time_stamp, msg);
			break;
		case QtFatalMsg:
			msg_level = 0;
			log_statement = QString("%1 - [Fatal] %2").arg(time_stamp, msg);
			break;
		case QtInfoMsg:
			msg_level = 1;
			log_statement = QString("%1 - [Info] %2").arg(time_stamp, msg);
			break;
		case QtWarningMsg:
			msg_level = 2;
			log_statement = QString("%1 - [Warning] %2").arg(time_stamp, msg);
			break;
		case QtDebugMsg:
		default:
			msg_level = 3;
			log_statement = QString("%1 - [Debug] %2").arg(time_stamp, msg);
	}

	// Log levels:
	// 0: only critical and fatal
	// 1: += info
	// 2: += warning
	// 3: += debug
	if (msg_level <= log_level)
	{
		printf("%s", qUtf8Printable(log_statement.replace("\"", "")));
		printf("\n");

		QTextStream out_stream(&gsvar_server_log_file);
		out_stream.setCodec("UTF-8");
		out_stream.setGenerateByteOrderMark(false);
		out_stream << log_statement << endl;
	}

	if (type == QtFatalMsg)
	{
		abort();
	}
}

int main(int argc, char **argv)
{
	gsvar_server_log_file.open(QIODevice::WriteOnly | QIODevice::Append);

	QCoreApplication app(argc, argv);

	QCommandLineParser parser;
	parser.setApplicationDescription("GSvar file server");
	parser.addHelpOption();
	parser.addVersionOption();
	QCommandLineOption serverPortOption(QStringList() << "p" << "port",
			QCoreApplication::translate("main", "Server port number"),
			QCoreApplication::translate("main", "port"));
	parser.addOption(serverPortOption);
	QCommandLineOption logLevelOption(QStringList() << "l" << "log",
			QCoreApplication::translate("main", "Log level"),
			QCoreApplication::translate("main", "logging"));
	parser.addOption(logLevelOption);
	parser.process(app);
	QString port = parser.value(serverPortOption);
	QString log_level_option = parser.value(logLevelOption);

	if (!log_level_option.isEmpty())
	{
		qInfo().noquote() << "Log level parameter has been provided through the command line arguments:" + log_level_option;
		log_level = log_level_option.toInt();
	}
	else {
		qInfo().noquote() << "Using log level from the application settings:" + QString::number(ServerHelper::getNumSettingsValue("log_level"));
		log_level = ServerHelper::getNumSettingsValue("log_level");
	}

	qInstallMessageHandler(interceptLogMessage);

	EndpointManager::appendEndpoint(Endpoint{
						"",
						QMap<QString, ParamProps>{},
						RequestMethod::GET,
						ContentType::TEXT_HTML,
						false,
						"Index page with general information",
						&EndpointHandler::serveIndexPage
					});
	EndpointManager::appendEndpoint(Endpoint{
						"favicon.ico",
						QMap<QString, ParamProps>{},
						RequestMethod::GET,
						ContentType::IMAGE_PNG,
						false,
						"Favicon to avoid warnings from the browser",
						&EndpointHandler::serveFavicon
					});
	EndpointManager::appendEndpoint(Endpoint{
						"info",
						QMap<QString, ParamProps>{},
						RequestMethod::GET,
						ContentType::APPLICATION_JSON,
						false,
						"General information about this API",
						&EndpointHandler::serveApiInfo
					});
	EndpointManager::appendEndpoint(Endpoint{
						"static",
						QMap<QString, ParamProps>{
						   {"filename", ParamProps{ParamProps::ParamCategory::PATH_PARAM, true, "Name of the file to be served"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_HTML,
						false,
						"Static content served from the server root folder (defined in the config file)",
						&EndpointController::serveStaticFromServerRoot
				   });

	EndpointManager::appendEndpoint(Endpoint{
						"static",
						QMap<QString, ParamProps>{
						   {"filename", ParamProps{ParamProps::ParamCategory::PATH_PARAM, false, "Name of the file to be served"}}
						},
						RequestMethod::HEAD,
						ContentType::TEXT_HTML,
						false,
						"Size of the static content served from the server root folder (defined in the config file)",
						&EndpointController::serveStaticFromServerRoot
				   });

	EndpointManager::appendEndpoint(Endpoint{
						"protected",
						QMap<QString, ParamProps>{
						   {"filename", ParamProps{ParamProps::ParamCategory::PATH_PARAM, true, "Name of the file to be served"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_HTML,
						true,
						"Protected static files",
						&EndpointController::serveStaticFromServerRoot
				   });

	EndpointManager::appendEndpoint(Endpoint{
						"cache",
						QMap<QString, ParamProps>{
						   {"filename", ParamProps{ParamProps::ParamCategory::PATH_PARAM, false, "Name of the file to be served"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_HTML,
						false,
						"Static content served from the server cache",
						&EndpointController::serveStaticFileFromCache
				   });
	EndpointManager::appendEndpoint(Endpoint{
						"temp",
						QMap<QString, ParamProps>{
						   {"id", ParamProps{ParamProps::ParamCategory::PATH_PARAM, false, "Unique id pointing to a file"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_HTML,
						false,
						"Static file served via secure temporary URL",
						&EndpointController::serveStaticForTempUrl
				   });

	EndpointManager::appendEndpoint(Endpoint{
						"temp",
						QMap<QString, ParamProps>{
							{"id", ParamProps{ParamProps::ParamCategory::PATH_PARAM, false, "Unique id pointing to a folder"}},
							{"filename", ParamProps{ParamProps::ParamCategory::PATH_PARAM, true, "Filename in a folder with a temporary URL"}}
						},
						RequestMethod::HEAD,
						ContentType::TEXT_HTML,
						false,
						"Size of the static file served via secure temporary URL",
						&EndpointController::serveStaticForTempUrl
				   });

	EndpointManager::appendEndpoint(Endpoint{
						"help",
						QMap<QString, ParamProps>{
						   {"endpoint", ParamProps{ParamProps::ParamCategory::PATH_PARAM, true,
							"Endpoint path the help is requested for. Help for all endpoints wiil be provided, if this parameter is ommited"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_HTML,
						false,
						"Help page on the usage of the endpoints",
						&EndpointController::serveEndpointHelp
					});

	EndpointManager::appendEndpoint(Endpoint{
						"file_location",
						QMap<QString, ParamProps> {
						   {"ps_url_id", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "An id of a temporary URL pointing to a specific processed sample"}},
						   {"type", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, true, "Format of the requested file(s)"}},
						   {"path", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, true, "Returns an absolute path on the server, if set to 'absolute'"}},
						   {"locus", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, true, "Locus for repeat expansion image"}},
						   {"return_if_missing", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, true, "Return file info, if the file is missing"}}
						},
						RequestMethod::GET,
						ContentType::APPLICATION_JSON,
						false,
						"Retrieve file location information for specific file types",
						&EndpointHandler::locateFileByType
					});

	EndpointManager::appendEndpoint(Endpoint{
						"project_file",
						QMap<QString, ParamProps> {
						   {"ps_id", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Processed sample id"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_PLAIN,
						false,
						"Temporary URL leading to a specific project file (GSvar file)",
						&EndpointHandler::locateProjectFile
					});

	EndpointManager::appendEndpoint(Endpoint{
						"project_file",
						QMap<QString, ParamProps> {
						   {"ps_url_id", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "An id of a temporary URL pointing to a specific processed sample"}}
						},
						RequestMethod::PUT,
						ContentType::APPLICATION_JSON,
						false,
						"Update an existing project file (GSvar file)",
						&EndpointHandler::saveProjectFile
					});

	EndpointManager::appendEndpoint(Endpoint{
						"file_info",
						QMap<QString, ParamProps> {
						   {"file", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Filename with its full path"}}
						},
						RequestMethod::GET,
						ContentType::APPLICATION_JSON,
						false,
						"Detailed information about a specific file",
						&EndpointController::getFileInfo
					});


	EndpointManager::appendEndpoint(Endpoint{
						"ps_regions",
						QMap<QString, ParamProps> {
						   {"sys_id", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Processing system id"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_PLAIN,
						false,
						"Processing system regions",
						&EndpointHandler::getProcessingSystemRegions
					});

	EndpointManager::appendEndpoint(Endpoint{
						"ps_amplicons",
						QMap<QString, ParamProps> {
						   {"sys_id", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Processing system id"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_PLAIN,
						false,
						"Processing system amplicons",
						&EndpointHandler::getProcessingSystemAmplicons
					});

	EndpointManager::appendEndpoint(Endpoint{
						"ps_genes",
						QMap<QString, ParamProps> {
						   {"sys_id", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Processing system id"}}
						},
						RequestMethod::GET,
						ContentType::TEXT_PLAIN,
						false,
						"Processing system genes",
						&EndpointHandler::getProcessingSystemGenes
					});

	EndpointManager::appendEndpoint(Endpoint{
						"secondary_analyses",
						QMap<QString, ParamProps> {
							{"ps_name", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Processed sample name"}},
							{"type", ParamProps{ParamProps::ParamCategory::GET_URL_PARAM, false, "Analysis type"}}
						},
						RequestMethod::GET,
						ContentType::APPLICATION_JSON,
						false,
						"Secondary analyses list",
						&EndpointHandler::getSecondaryAnalyses
					});

	int port_number = ServerHelper::getNumSettingsValue("server_port");

	if (!port.isEmpty())
	{
		qInfo() << "Server port has been provided through the command line arguments:" + port;
		port_number = port.toInt();
	}
	else {
		qInfo() << "Using port number from the application settings";
	}

	qInfo() << "SSL version used for build: " << QSslSocket::sslLibraryBuildVersionString();

	HttpsServer sslserver(port_number);
	return app.exec();
}
