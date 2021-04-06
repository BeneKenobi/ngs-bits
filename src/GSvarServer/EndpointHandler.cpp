#include "EndpointHandler.h"

EndpointHandler::EndpointHandler()
{
}

bool EndpointHandler::isValidUser(QString name, QString password)
{
	try
	{
		NGSD db;
		QString message = db.checkPassword(name, password, true);
		if (message.isEmpty())
		{
			return true;
		}
		else
		{
			return false;
		}

	}
	catch (DatabaseException& e)
	{
		qCritical() << e.message();
	}
	return false;
}

QList<QString> EndpointHandler::getAnalysisFiles(QString sample_name, bool search_multi)
{
	QList<QString> files;
	try
	{
		//convert name to file
		NGSD db;
		QString processed_sample_id = db.processedSampleId(sample_name);
		QString analysis_file = db.processedSamplePath(processed_sample_id, PathType::GSVAR);

		//determine all analyses of the sample
		if (QFile::exists(analysis_file)) files << analysis_file;

		//somatic tumor sample > ask user if he wants to open the tumor-normal pair
		QString normal_sample = db.normalSample(processed_sample_id);
		if (normal_sample!="")
		{
			files << db.secondaryAnalyses(sample_name + "-" + normal_sample, "somatic");
		}
		//check for germline trio/multi analyses
		else if (search_multi)
		{
			files << db.secondaryAnalyses(sample_name, "trio");
			files << db.secondaryAnalyses(sample_name, "multi sample");
		}
	}
	catch (Exception& e)
	{
		qWarning() << "Error opening processed sample from NGSD:" << e.message();
	}
	return files;
}

HttpResponse EndpointHandler::serveIndexPage(HttpRequest request)
{	
	return EndpointHelper::serveStaticFile(":/assets/client/info.html", ByteRange{}, ContentType::TEXT_HTML, false);
}

HttpResponse EndpointHandler::serveApiInfo(HttpRequest request)
{
	return EndpointHelper::serveStaticFile(":/assets/client/api.json", ByteRange{}, ContentType::APPLICATION_JSON, false);
}

HttpResponse EndpointHandler::serveTempUrl(HttpRequest request)
{
	UrlEntity url_entity = UrlManager::getURLById(request.getPathParams()[0]);
	if (url_entity.filename_with_path.isEmpty())
	{
		return HttpResponse(HttpError{StatusCode::NOT_FOUND, request.getContentType(), "Could not find a file linked to the id: " + request.getPathParams()[0]});
	}

	qDebug() << "Serving file: " << url_entity.filename_with_path;
//	return EndpointHelper::serveStaticFile(url_entity.filename_with_path, ByteRange{}, HttpProcessor::getContentTypeByFilename(url_entity.filename_with_path), false);
	return EndpointHelper::streamStaticFile(url_entity.filename_with_path, false);
}

HttpResponse EndpointHandler::locateFileByType(HttpRequest request)
{
	qDebug() << "File location service";
	if (!request.getUrlParams().contains("ps"))
	{
		return HttpResponse(HttpError{StatusCode::BAD_REQUEST, request.getContentType(), "Sample id has not been provided"});
	}
	QString ps = request.getUrlParams().value("ps");
	qDebug() << "PS" << ps;
	QString found_file;
	if (ps.indexOf("/")>-1){

		UrlEntity url_entity = UrlManager::getURLById(ps.split("/").last().trimmed());
		qDebug() << "FOUND PROJECT" << url_entity.filename_with_path;
		found_file = url_entity.filename_with_path;
	}
//	else
//	{
//		found_file = getGSvarFile(request.getUrlParams().value("ps"), false);
//	}

	bool return_if_missing = true;
	if (!request.getUrlParams().contains("return_if_missing"))
	{
		if (request.getUrlParams().value("return_if_missing") == "0")
		{
			return_if_missing = false;
		}
	}

	if (found_file.isEmpty())
	{
		return HttpResponse(HttpError{StatusCode::NOT_FOUND, request.getContentType(), "Could not find the sample: " + request.getUrlParams().value("ps")});
	}

	VariantList variants;
	variants.load(found_file);

	FileLocationProviderLocal* file_locator = new FileLocationProviderLocal(found_file, variants.getSampleHeader(), variants.type());

	qDebug() << "found_file " << found_file;
	QList<FileLocation> file_list {};
	QJsonDocument json_doc_output {};
	QJsonArray json_list_output {};

	if(request.getUrlParams()["type"].toLower() == "analysisvcf")
	{
		file_list << file_locator->getAnalysisVcf();
	}
	else if(request.getUrlParams()["type"].toLower() == "analysissv")
	{
		file_list << file_locator->getAnalysisSvFile();
	}
	else if(request.getUrlParams()["type"].toLower() == "analysiscnv")
	{
		file_list << file_locator->getAnalysisCnvFile();
	}
	else if(request.getUrlParams()["type"].toLower() == "analysismosaiccnv")
	{
		file_list << file_locator->getAnalysisMosaicCnvFile();
	}
	else if(request.getUrlParams()["type"].toLower() == "analysisupd")
	{
		file_list << file_locator->getAnalysisUpdFile();
	}
	else if(request.getUrlParams()["type"].toLower() == "repeatexpansionimage")
	{
		if (!request.getUrlParams().contains("locus"))
		{
			return HttpResponse(HttpError{StatusCode::BAD_REQUEST, request.getContentType(), "Locus value has not been provided"});
		}
		file_list << file_locator->getRepeatExpansionImage(request.getUrlParams().value("locus"));
	}
	else if(request.getUrlParams()["type"].toLower() == "bam")
	{
		file_list = file_locator->getBamFiles(return_if_missing);
	}
	else if(request.getUrlParams()["type"].toLower() == "cnvcoverage")
	{
		file_list = file_locator->getCnvCoverageFiles(return_if_missing);
	}
	else if(request.getUrlParams()["type"].toLower() == "baf")
	{
		file_list = file_locator->getBafFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "mantaevidence")
	{
		file_list = file_locator->getMantaEvidenceFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "circosplot")
	{
		file_list = file_locator->getCircosPlotFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "vcf")
	{
		file_list = file_locator->getVcfFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "repeatexpansion")
	{
		file_list = file_locator->getRepeatExpansionFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "prs")
	{
		file_list = file_locator->getPrsFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "lowcoverage")
	{
		file_list = file_locator->getLowCoverageFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "copynumbercall")
	{
		file_list = file_locator->getCopyNumberCallFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "roh")
	{
		file_list = file_locator->getRohFiles(return_if_missing);
	}
	else if (request.getUrlParams()["type"].toLower() == "somaticcnvcoverage")
	{
		file_list << file_locator->getSomaticCnvCoverageFile();
	}
	else if (request.getUrlParams()["type"].toLower() == "somaticcnvcall")
	{
		file_list << file_locator->getSomaticCnvCallFile();
	}
	else if (request.getUrlParams()["type"].toLower() == "somaticlowcoverage")
	{
		file_list << file_locator->getSomaticLowCoverageFile();
	}
	else if (request.getUrlParams()["type"].toLower() == "somaticmsi")
	{
		file_list << file_locator->getSomaticMsiFile();
	}
	else
	{
		FileLocation gsvar_file(
			request.getUrlParams()["ps"],
			PathType::GSVAR,
			found_file,
			true
		);
		file_list.append(gsvar_file);
	}

	for (int i = 0; i < file_list.count(); ++i)
	{
		qDebug() << file_list[i].filename;
		QJsonObject cur_json_item {};
		cur_json_item.insert("id", file_list[i].id);
		cur_json_item.insert("type", FileLocation::typeToString(file_list[i].type));
		bool needs_url = true;
		if (request.getUrlParams().contains("path"))
		{
			if (request.getUrlParams()["path"].toLower() == "absolute") needs_url = false;

		}
		if (needs_url)
		{
			cur_json_item.insert("filename", createTempUrl(file_list[i].filename));
		}
		else
		{
			cur_json_item.insert("filename", file_list[i].filename);
		}

		cur_json_item.insert("exists", file_list[i].exists);
		json_list_output.append(cur_json_item);
	}

	json_doc_output.setArray(json_list_output);
	return HttpResponse(false, false, "", EndpointHelper::generateHeaders(json_doc_output.toJson().length(), ContentType::APPLICATION_JSON), json_doc_output.toJson());
}

HttpResponse EndpointHandler::locateProjectFile(HttpRequest request)
{
	qDebug() << "Project file location";
	QJsonDocument json_doc_output;
	QJsonArray json_list_output;

	bool search_multi = false;
	if (request.getUrlParams().contains("multi"))
	{
		if (request.getUrlParams().value("multi") == "1")
		{
			search_multi = true;
		}
	}

	QList<QString> found_files = getAnalysisFiles(request.getUrlParams().value("ps"), search_multi);
	for (int i = 0; i < found_files.count(); i++)
	{
		json_list_output.append(createTempUrl(found_files[i]));
	}
	json_doc_output.setArray(json_list_output);
	return HttpResponse(false, false, "", EndpointHelper::generateHeaders(json_doc_output.toJson().length(), ContentType::APPLICATION_JSON), json_doc_output.toJson());
}

HttpResponse EndpointHandler::performLogin(HttpRequest request)
{
	QByteArray body {};
	if (!request.getFormUrlEncoded().contains("name") || !request.getFormUrlEncoded().contains("password"))
	{
		return HttpResponse(HttpError{StatusCode::FORBIDDEN, request.getContentType(), "No username or/and password were found"});
	}

	if (isValidUser(request.getFormUrlEncoded()["name"], request.getFormUrlEncoded()["password"]))
	{
		QString secure_token = ServerHelper::generateUniqueStr();
		Session cur_session = Session{request.getFormUrlEncoded()["name"], QDateTime::currentDateTime()};

		SessionManager::addNewSession(secure_token, cur_session);
		body = secure_token.toLocal8Bit();
		return HttpResponse{false, false, "", EndpointHelper::generateHeaders(body.length(), ContentType::TEXT_PLAIN), body};
	}

	return HttpResponse(HttpError{StatusCode::UNAUTHORIZED, request.getContentType(), "Invalid username or password"});
}

HttpResponse EndpointHandler::performLogout(HttpRequest request)
{
	QByteArray body {};
	if (!request.getFormUrlEncoded().contains("token"))
	{
		return HttpResponse(HttpError{StatusCode::FORBIDDEN, request.getContentType(), "Secure token has not been provided"});
	}
	if (SessionManager::isTokenValid(request.getFormUrlEncoded()["token"]))
	{
		try
		{
			SessionManager::removeSession(request.getFormUrlEncoded()["token"]);
		} catch (Exception& e)
		{
			return HttpResponse(HttpError{StatusCode::INTERNAL_SERVER_ERROR, request.getContentType(), e.message()});
		}
		body = request.getFormUrlEncoded()["token"].toLocal8Bit();
		return HttpResponse{false, false, "", EndpointHelper::generateHeaders(body.length(), ContentType::TEXT_PLAIN), body};
	}
	return HttpResponse(HttpError{StatusCode::FORBIDDEN, request.getContentType(), "You have provided an invalid token"});
}

QString EndpointHandler::createTempUrl(QString filename)
{
	QString id = ServerHelper::generateUniqueStr();
	UrlManager::addUrlToStorage(id, filename);
	return ServerHelper::getStringSettingsValue("server_host") +
			+ ":" + ServerHelper::getStringSettingsValue("server_port") +
			+ "/v1/temp/" + id;
}
