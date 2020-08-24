#include "TestFramework.h"
#include "Settings.h"
#include "TestFrameworkNGS.h"

TEST_CLASS(VcfLeftNormalize_Test)
{
Q_OBJECT
private slots:
	
	void test_01()
	{
		QString ref_file = Settings::string("reference_genome", true);
		if (ref_file=="") SKIP("Test needs the reference genome!");

		EXECUTE("VcfLeftNormalize", "-in " + TESTDATA("data_in/VcfLeftNormalize_in1.vcf") + " -out out/VcfLeftNormalize_out1.vcf -comp 0 -ref " + ref_file);
		COMPARE_FILES("out/VcfLeftNormalize_out1.vcf", TESTDATA("data_out/VcfLeftNormalize_out1.vcf"));
		VCF_IS_VALID("out/VcfLeftNormalize_out1.vcf")
	}

};

