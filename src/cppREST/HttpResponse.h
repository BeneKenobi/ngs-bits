#ifndef HTTPRESPONSE_H
#define HTTPRESPONSE_H

#include "cppREST_global.h"
#include <QByteArray>
#include <QString>
#include <QDir>
#include <QJsonDocument>
#include <QJsonObject>
#include "HttpParts.h"
#include "HttpProcessor.h"
#include "HtmlEngine.h"

class CPPRESTSHARED_EXPORT HttpResponse : public QByteArray
{
public:
	HttpResponse();
	HttpResponse(QByteArray response_data);
	HttpResponse(BasicResponseData data);
	HttpResponse(BasicResponseData data, QByteArray payload);
	HttpResponse(ResponseStatus status, ContentType content_type, QString message);

	void setIsStream(bool is_stream);
	bool isStream();

	void setFilename(QString filename);
	QString getFilename();

	void setStatus(ResponseStatus response_status);
	ResponseStatus getStatus();

	QByteArray getStatusLine();

	int getStatusCode();

	void setHeaders(QByteArray headers);
	void addHeader(QString header);
	QByteArray getHeaders();

	void setPayload(QByteArray payload);
	QByteArray getPayload();

	void setRangeNotSatisfiableHeaders(BasicResponseData data);

private:
	void readBasicResponseData(BasicResponseData data);
	QByteArray generateRegularHeaders(BasicResponseData data);
	QByteArray generateChunkedStreamHeaders(BasicResponseData data);
	QByteArray generateRangeNotSatisfiableHeaders(BasicResponseData data);
	QString getFileNameWithExtension(QString filename_with_path);

protected:
	bool is_stream_;
	QString filename_;
	ResponseStatus response_status_;
	QByteArray status_line_;
	QByteArray headers_;
	QByteArray payload_;
	int getContentLength();
	void updateResponseData();
};

#endif // HTTPRESPONSE_H
