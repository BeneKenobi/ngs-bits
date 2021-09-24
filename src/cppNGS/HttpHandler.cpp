#include "HttpHandler.h"
#include "Exceptions.h"
#include "Settings.h"
#include "Helper.h"

#include <QNetworkRequest>
#include <QNetworkReply>
#include <QEventLoop>
#include <QNetworkProxy>
#include <QAuthenticator>
#include <QFile>
#include <QPointer>
#include <QHttpMultiPart>

HttpHandler::HttpHandler(HttpRequestHandler::ProxyType proxy_type, QObject* parent)
	: QObject(parent)
	, nmgr_()
	, proxy_type_(proxy_type)
	, headers_()
{
	//default headers
	setHeader("User-Agent", "GSvar");
	setHeader("X-Custom-User-Agent", "GSvar");
}

const HttpHeaders& HttpHandler::headers() const
{
	return headers_;
}

void HttpHandler::setHeader(const QByteArray& key, const QByteArray& value)
{
	headers_.insert(key, value);
}

QByteArray HttpHandler::get(QString url, const HttpHeaders& add_headers)
{
	return HttpRequestHandler(proxy_type_, 0).get(url, add_headers);
}

QByteArray HttpHandler::post(QString url, const QByteArray& data, const HttpHeaders& add_headers)
{
	return HttpRequestHandler(proxy_type_, 0).post(url, data, add_headers);
}

QByteArray HttpHandler::post(QString url, QHttpMultiPart* parts, const HttpHeaders& add_headers)
{
	return HttpRequestHandler(proxy_type_, 0).post(url, parts, add_headers);
}



