#include "TestFramework.h"
#include "Settings.h"
#include "NGSD.h"

TEST_CLASS(NGSDImportORPHA_Test)
{
Q_OBJECT
private slots:
	
	void default_parameters()
	{
		QString host = Settings::string("ngsd_test_host");
		if (host=="") SKIP("Test needs access to the NGSD test database!");

		//init
		NGSD db(true);
		db.init();

		//test
		EXECUTE("NGSDImportORPHA", "-test -terms " + TESTDATA("data_in/NGSDImportORPHA_terms.xml") + " -genes " + TESTDATA("data_in/NGSDImportORPHA_genes.xml"));

		//check
		I_EQUAL(db.getValue("SELECT count(*) FROM disease_term").toInt(), 2)
		I_EQUAL(db.getValue("SELECT count(*) FROM disease_gene").toInt(), 3)
	}
};

