#ifndef ENDPOINTHELPER_H
#define ENDPOINTHELPER_H

#include "cppREST_global.h"
#include "HttpParts.h"
#include "HtmlEngine.h"
#include "SessionManager.h"
#include "FileCache.h"
#include "EndpointManager.h"
#include "UrlManager.h"
#include "HttpResponse.h"


class CPPRESTSHARED_EXPORT EndpointHelper
{
public:
	static bool isEligibileToAccess(HttpRequest request);
	static QString getFileNameWithExtension(QString filename_with_path);
	static StaticFile readFileContent(QString filename, ByteRange byte_range);
	static QString addFileToCache(QString filename);
	static HttpResponse serveStaticFile(QString filename, ByteRange byte_range, ContentType type, bool is_downloadable);
	static HttpResponse serveStaticFileFromCache(QString id, ByteRange byte_range, ContentType type, bool is_downloadable);
	static HttpResponse streamStaticFile(QString filename, bool is_downloadable);
	static HttpResponse serveFolderContent(QString folder);
	static QByteArray generateHeaders(QString filename, int length, ByteRange byte_range, qint64 file_size, ContentType type, bool is_downloadable);
	static QByteArray generateHeaders(int length, ContentType type);
	static HttpResponse listFolderContent(HttpRequest request);
	static HttpResponse serveEndpointHelp(HttpRequest request);
	static HttpResponse serveStaticFile(HttpRequest request);
	static HttpResponse serveStaticFileFromCache(HttpRequest request);
	static HttpResponse streamStaticFile(HttpRequest request);
	static HttpResponse serveProtectedStaticFile(HttpRequest request);
	static HttpResponse getFileInfo(HttpRequest request);
	static HttpResponse serveFolderListing(QList<FolderItem> items);

protected:
	EndpointHelper();	

private:
	static EndpointHelper& instance();
};

#endif // ENDPOINTHELPER_H
