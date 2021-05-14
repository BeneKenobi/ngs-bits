#include "RequestParser.h"

RequestPaser::RequestPaser(QByteArray *request, QString client_address)
	: raw_request_(request)
	, client_address_(client_address)
{
}

HttpRequest RequestPaser::getRequest()
{
	HttpRequest request;
	QList<QByteArray> body = splitRequestIntoLines();

	request.setRemoteAddress(client_address_);

	for (int i = 0; i < body.count(); ++i)
	{
		// First line with method type and URL
		if (i == 0)
		{
			QList<QByteArray> request_info = body[i].split(' ');
			if (request_info.length() < 2)
			{
				THROW(Exception, "Cannot process the request. It is possible the URL is missing or incorrect");
				return request;
			}
			request.setMethod(inferRequestMethod(request_info[0].toUpper()));

			QList<QString> path_items = QString(request_info[1]).split('/');
			qDebug() << "path_items = " << path_items;
			qDebug() << "";
			request.setPrefix(getRequestPrefix(path_items));
			request.setPath(getRequestPath(path_items));
			request.setPathParams(getRequestPathParams(path_items));
			request.setUrlParams(getVariables(getVariableSequence(request_info[1])));
			continue;
		}

		// Reading headers and params
		int header_separator = body[i].indexOf(':');
		int param_separator = body[i].indexOf('=');
		if ((header_separator == -1) && (param_separator == -1) && (body[i].length() > 0))
		{
			THROW(Exception, "Malformed element: " + body[i]);
			return request;
		}

		if (header_separator > -1)
		{
			request.addHeader(body[i].left(header_separator).toLower(), body[i].mid(header_separator+1).trimmed());
		}
		else if (param_separator > -1)
		{
			request.setFormUrlEncoded(getVariables(body[i]));
		}
	}

	request.setContentType(ContentType::TEXT_HTML);
	if (request.getHeaders().contains("accept"))
	{
		if (HttpProcessor::getContentTypeFromString(request.getHeaders()["accept"]) == ContentType::APPLICATION_JSON)
		{
			request.setContentType(ContentType::APPLICATION_JSON);
		}
	}

	return request;
}

QList<QByteArray> RequestPaser::splitRequestIntoLines()
{
	QList<QByteArray> request_items;
	request_items = raw_request_->split('\n');
	for (int i = 0; i < request_items.count(); ++i)
	{
		request_items[i] = request_items[i].trimmed();
	}
	return request_items;
}

QList<QByteArray> RequestPaser::getKeyValuePair(QByteArray input)
{
	QList<QByteArray> result;

	if (input.indexOf('=')>-1)
	{
		result = input.split('=');
	}

	return result;
}

QMap<QString, QString> RequestPaser::getVariables(QByteArray input)
{
	QMap<QString, QString> url_vars {};
	QList<QByteArray> var_list = input.split('&');
	QByteArray cur_key {};

	for (int i = 0; i < var_list.count(); ++i)
	{
		QList<QByteArray> pair = getKeyValuePair(var_list[i]);
		if (pair.length()==2)
		{
			url_vars.insert(pair[0], pair[1]);
		}
	}

	return url_vars;
}

QByteArray RequestPaser::getVariableSequence(QByteArray url)
{
	QByteArray var_string;
	if (url.indexOf('?') == -1) return var_string;
	return url.split('?')[1];
}

QString RequestPaser::getRequestPrefix(QList<QString> path_items)
{
	if (path_items.count()>1)
	{
		return ServerHelper::getUrlWithoutParams(path_items[1]);
	}
	return "";
}

QString RequestPaser::getRequestPath(QList<QString> path_items)
{
	if (path_items.count()>2)
	{
		return ServerHelper::getUrlWithoutParams(path_items[2]);
	}
	return "";
}

QList<QString> RequestPaser::getRequestPathParams(QList<QString> path_items)
{
	QList<QString> params {};
	if (path_items.count()>3)
	{
		for (int p = 3; p < path_items.count(); ++p)
		{
			if (!path_items[p].trimmed().isEmpty())
			{
				params.append(path_items[p].trimmed());
			}
		}
	}
	return params;
}

RequestMethod RequestPaser::inferRequestMethod(QByteArray input)
{
	if (input.toUpper() == QByteArrayLiteral("GET"))
	{
		return RequestMethod::GET;
	}
	if (input.toUpper() == QByteArrayLiteral("POST"))
	{
		return RequestMethod::POST;
	}
	if (input.toUpper() == QByteArrayLiteral("DELETE"))
	{
		return RequestMethod::DELETE;
	}
	if (input.toUpper() == QByteArrayLiteral("PUT"))
	{
		return RequestMethod::PUT;
	}
	if (input.toUpper() == QByteArrayLiteral("PATCH"))
	{
		return RequestMethod::PATCH;
	}
	if (input.toUpper() == QByteArrayLiteral("HEAD"))
	{
		return RequestMethod::HEAD;
	}

	THROW(ArgumentException, "Incorrect request method");
}
