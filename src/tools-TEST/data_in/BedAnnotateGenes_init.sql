
INSERT INTO `gene` (`id`, `hgnc_id`, `symbol`, `name`, `type`) VALUES
(81676, 1101, 'BRCA2', 'breast cancer 2, early onset', 'protein-coding gene'),
(103240, 9121, 'PMS1', 'PMS1 homolog 1, mismatch repair system component', 'protein-coding gene');

INSERT INTO `gene_transcript` (`id`, `gene_id`, `name`, `version`, `source`, `chromosome`, `start_coding`, `end_coding`, `strand`) VALUES
(377820, 81676, 'CCDS9344', '1', 'ccds', '13', 32890598, 32972907, '+'),
(377819, 81676, 'ENST00000380152', '1', 'ensembl', '13', 32890598, 32972907, '+'),
(377821, 81676, 'ENST00000544455', '1', 'ensembl', '13', 32890598, 32972907, '+'),
(409247, 103240, 'CCDS2302', '1', 'ccds', '2', 190656536, 190742162, '+'),
(409249, 103240, 'CCDS46473', '1', 'ccds', '2', 190656536, 190742162, '+'),
(409242, 103240, 'CCDS46474', '1', 'ccds', '2', 190656536, 190742162, '+'),
(409240, 103240, 'ENST00000374826', '1', 'ensembl', '2', 190656536, 190671228, '+'),
(409241, 103240, 'ENST00000409823', '1', 'ensembl', '2', 190656536, 190742162, '+'),
(409243, 103240, 'ENST00000409985', '1', 'ensembl', '2', 190656536, 190670608, '+'),
(409244, 103240, 'ENST00000418224', '1', 'ensembl', '2', 190682853, 190742162, '+'),
(409245, 103240, 'ENST00000432292', '1', 'ensembl', '2', 190682853, 190742162, '+'),
(409246, 103240, 'ENST00000441310', '1', 'ensembl', '2', 190656536, 190742162, '+'),
(409248, 103240, 'ENST00000447232', '1', 'ensembl', '2', 190656536, 190742162, '+');

INSERT INTO `gene_exon` (`transcript_id`, `start`, `end`) VALUES
(377819, 32889611, 32889804),
(377819, 32890559, 32890664),
(377819, 32893214, 32893462),
(377819, 32899213, 32899321),
(377819, 32900238, 32900287),
(377819, 32900379, 32900419),
(377819, 32900636, 32900750),
(377819, 32903580, 32903629),
(377819, 32905056, 32905167),
(377819, 32906409, 32907524),
(377819, 32910402, 32915333),
(377819, 32918695, 32918790),
(377819, 32920964, 32921033),
(377819, 32928998, 32929425),
(377819, 32930565, 32930746),
(377819, 32931879, 32932066),
(377819, 32936660, 32936830),
(377819, 32937316, 32937670),
(377819, 32944539, 32944694),
(377819, 32945093, 32945237),
(377819, 32950807, 32950928),
(377819, 32953454, 32953652),
(377819, 32953887, 32954050),
(377819, 32954144, 32954282),
(377819, 32968826, 32969070),
(377819, 32971035, 32971181),
(377819, 32972299, 32973347),
(377820, 32890598, 32890664),
(377820, 32893214, 32893462),
(377820, 32899213, 32899321),
(377820, 32900238, 32900287),
(377820, 32900379, 32900419),
(377820, 32900636, 32900750),
(377820, 32903580, 32903629),
(377820, 32905056, 32905167),
(377820, 32906409, 32907524),
(377820, 32910402, 32915333),
(377820, 32918695, 32918790),
(377820, 32920964, 32921033),
(377820, 32928998, 32929425),
(377820, 32930565, 32930746),
(377820, 32931879, 32932066),
(377820, 32936660, 32936830),
(377820, 32937316, 32937670),
(377820, 32944539, 32944694),
(377820, 32945093, 32945237),
(377820, 32950807, 32950928),
(377820, 32953454, 32953652),
(377820, 32953887, 32954050),
(377820, 32954144, 32954282),
(377820, 32968826, 32969070),
(377820, 32971035, 32971181),
(377820, 32972299, 32972907),
(377821, 32889617, 32889804),
(377821, 32890559, 32890664),
(377821, 32893214, 32893462),
(377821, 32899213, 32899321),
(377821, 32900238, 32900287),
(377821, 32900379, 32900419),
(377821, 32900636, 32900750),
(377821, 32903580, 32903629),
(377821, 32905056, 32905167),
(377821, 32906409, 32907524),
(377821, 32910402, 32915333),
(377821, 32918695, 32918790),
(377821, 32920964, 32921033),
(377821, 32928998, 32929425),
(377821, 32930565, 32930746),
(377821, 32931879, 32932066),
(377821, 32936660, 32936830),
(377821, 32937316, 32937670),
(377821, 32944539, 32944694),
(377821, 32945093, 32945237),
(377821, 32950807, 32950928),
(377821, 32953454, 32953652),
(377821, 32953887, 32954050),
(377821, 32954144, 32954282),
(377821, 32968826, 32969070),
(377821, 32971035, 32971181),
(377821, 32972299, 32973347),
(377821, 32973746, 32973805),
(409240, 190649234, 190649319),
(409240, 190656516, 190656667),
(409240, 190660495, 190660677),
(409240, 190670378, 190670480),
(409240, 190671149, 190671572),
(409241, 190649227, 190649319),
(409241, 190656516, 190656667),
(409241, 190660495, 190660677),
(409241, 190670378, 190670480),
(409241, 190682743, 190682906),
(409241, 190717381, 190717503),
(409241, 190718665, 190718808),
(409241, 190718965, 190719854),
(409241, 190728469, 190728954),
(409241, 190732525, 190732655),
(409241, 190738222, 190738382),
(409241, 190741998, 190742162),
(409242, 190656536, 190656667),
(409242, 190660495, 190660677),
(409242, 190670378, 190670480),
(409242, 190682743, 190682906),
(409242, 190717381, 190717503),
(409242, 190718665, 190718808),
(409242, 190718965, 190719854),
(409242, 190728469, 190728954),
(409242, 190732525, 190732655),
(409242, 190738222, 190738382),
(409242, 190741998, 190742162),
(409243, 190649176, 190649319),
(409243, 190656516, 190656667),
(409243, 190660495, 190660677),
(409243, 190670378, 190671780),
(409244, 190649215, 190649319),
(409244, 190656516, 190656667),
(409244, 190660495, 190660677),
(409244, 190682743, 190682906),
(409244, 190708690, 190708806),
(409244, 190717381, 190717503),
(409244, 190718665, 190718808),
(409244, 190718965, 190719854),
(409244, 190728469, 190728954),
(409244, 190732525, 190732655),
(409244, 190738222, 190738382),
(409244, 190741998, 190742353),
(409245, 190649291, 190649319),
(409245, 190656516, 190656667),
(409245, 190682743, 190682906),
(409245, 190708690, 190708806),
(409245, 190717381, 190717503),
(409245, 190718665, 190718808),
(409245, 190718965, 190719854),
(409245, 190728469, 190728954),
(409245, 190732525, 190732655),
(409245, 190738222, 190738382),
(409245, 190741998, 190742259),
(409246, 190649107, 190649319),
(409246, 190656516, 190656667),
(409246, 190660495, 190660677),
(409246, 190670378, 190670480),
(409246, 190682743, 190682906),
(409246, 190708690, 190708806),
(409246, 190717381, 190717503),
(409246, 190718665, 190718808),
(409246, 190718965, 190719854),
(409246, 190728469, 190728954),
(409246, 190732525, 190732655),
(409246, 190738222, 190738382),
(409246, 190741998, 190742355),
(409247, 190656536, 190656667),
(409247, 190660495, 190660677),
(409247, 190670378, 190670480),
(409247, 190682743, 190682906),
(409247, 190708690, 190708806),
(409247, 190717381, 190717503),
(409247, 190718665, 190718808),
(409247, 190718965, 190719854),
(409247, 190728469, 190728954),
(409247, 190732525, 190732655),
(409247, 190738222, 190738382),
(409247, 190741998, 190742162),
(409248, 190649270, 190649319),
(409248, 190656516, 190656667),
(409248, 190660495, 190660677),
(409248, 190670378, 190670480),
(409248, 190682743, 190682906),
(409248, 190708690, 190708806),
(409248, 190717381, 190717503),
(409248, 190718665, 190718808),
(409248, 190718965, 190719854),
(409248, 190732525, 190732655),
(409248, 190738222, 190738382),
(409248, 190741998, 190742355),
(409249, 190656536, 190656667),
(409249, 190660495, 190660677),
(409249, 190670378, 190670480),
(409249, 190682743, 190682906),
(409249, 190708690, 190708806),
(409249, 190717381, 190717503),
(409249, 190718665, 190718808),
(409249, 190718965, 190719854),
(409249, 190732525, 190732655),
(409249, 190738222, 190738382),
(409249, 190741998, 190742162);
