#ifndef WEBSERVER_H
#define WEBSERVER_H

#include "cppREST_global.h"
#include <QObject>
#include <QDebug>
#include <QFile>
#include <QSslCertificate>
#include <QSslKey>
#include <QSslConfiguration>
#include <QSslSocket>
#include <QStandardPaths>
#include <QTimer>

#include "SslServer.h"
#include "UrlManager.h"

class CPPRESTSHARED_EXPORT WebServer : public QObject
{
    Q_OBJECT

public:
	WebServer(const quint16& port, const bool& insecure = false);

private:
	SslServer *server_;
};

#endif // WEBSERVER_H
