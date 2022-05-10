
-- device
INSERT INTO device (id, type, name) VALUES (1, 'HiSeq2500', 'Morpheus');

-- sequencing_run
INSERT INTO sequencing_run (id, name, fcid, device_id, recipe, quality) VALUES (1, 'First run', 'ABC', 1, '100+8+8+100', 'good');
INSERT INTO sequencing_run (id, name, fcid, device_id, recipe, quality) VALUES (2, 'Second run', 'XYZ', 1, '100+8+100', 'good');

-- user
INSERT INTO user (id, user_id, password, user_role, name, email, created, active) VALUES (99, 'ahuser', 's2d12kjg234hla0830t6hp9h3tt3t3tsdfg', 'user', 'The user', 'u@s.er', NOW(), '1');

-- sender
INSERT INTO sender (id, name) VALUES (1, 'sender');

-- project
INSERT INTO project (id, name, type, internal_coordinator_id, analysis) VALUES (1, 'First project', 'research', 1, 'variants');
INSERT INTO project (id, name, type, internal_coordinator_id, analysis) VALUES (2, 'Second project', 'diagnostic', 1, 'variants');

-- processing_system
INSERT INTO processing_system (id, name_manufacturer, shotgun, name_short, genome_id) VALUES (1, 'HaloPlex System', '1', 'hpSYSv1', 1);
INSERT INTO processing_system (id, name_manufacturer, shotgun, name_short, genome_id) VALUES 
(2, 'SureSelect Human All Exon v5', '1', 'ssHAEv5', 1);

-- sample
INSERT INTO sample (id, name, sample_type, species_id, gender, tumor, ffpe, sender_id, quality, disease_group, disease_status) VALUES
(1, 'NA12878', 'DNA', 1, 'female', '0', '0', 1, 'good', 'Diseases of the skin', 'Affected'),
(2, 'NA12879', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the skin', 'Affected'),
(3, 'NA12880', 'DNA', 1, 'female', '0', '0', 1, 'good', 'Diseases of the skin', 'Affected'),
(4, 'DUMMY', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(5, 'DUMMY2', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(6, 'DUMMY3', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(7, 'DUMMY4', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the immune system', 'Affected'),
(8, 'DUMMY5', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the immune system', 'Affected'),
(9, 'DUMMY6', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(10, 'DUMMY7', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(11, 'DUMMY8', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the visual system', 'Affected'),
(12, 'DUMMY9', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(13, 'DUMMY10', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Neoplasms', 'Affected'),
(14, 'DUMMY11', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the visual system', 'Affected'),
(15, 'DUMMY12', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the musculoskeletal system or connective tissue', 'Affected'),
(16, 'DUMMY13', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the musculoskeletal system or connective tissue', 'Affected'),
(17, 'DUMMY14', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the circulatory system', 'Affected'),
(18, 'DUMMY15', 'DNA', 1, 'male', '0', '0', 1, 'good', 'Diseases of the circulatory system', 'Affected');


-- processed_sample
INSERT INTO processed_sample (id, sample_id, process_id, sequencing_run_id, lane, operator_id, processing_system_id, project_id) VALUES
(1, 1, 1, 1, 1, 2, 1, 1),
(2, 2, 2, 2, 1, 2, 2, 2),
(3, 3, 3, 2, 1, 2, 2, 2),
(4, 4, 4, 2, 1, 2, 2, 2),
(5, 5, 5, 2, 1, 2, 2, 2),
(6, 6, 6, 2, 1, 2, 2, 2),
(7, 7, 7, 2, 1, 2, 2, 2),
(8, 8, 8, 2, 1, 2, 2, 2),
(9, 9, 9, 2, 1, 2, 2, 2),
(10, 10, 10, 2, 1, 2, 2, 2),
(11, 11, 11, 2, 1, 2, 2, 2),
(12, 12, 12, 2, 1, 2, 2, 2),
(13, 13, 13, 2, 1, 2, 2, 2),
(14, 14, 14, 2, 1, 2, 2, 2),
(15, 15, 15, 2, 1, 2, 2, 2),
(16, 16, 16, 2, 1, 2, 2, 2),
(17, 17, 17, 2, 1, 2, 2, 2),
(18, 18, 18, 2, 1, 2, 2, 2);
 
--variant
INSERT INTO variant (id, chr, start, end, ref, obs, comment, gnomad) VALUES
(1, 'chr1', 62247552, 62247552, 'C', 'G', 'from NA12880', 0.05),
(2, 'chr1', 62247574, 62247574, 'G', 'A', '', 0.04),
(3, 'chr1', 62263112, 62263112, 'A', 'G', '', 0.04),
(4, 'chr1', 62263166, 62263166, 'T', 'C', 'from DUMMY', 0.04),
(5, 'chr1', 119996708, 119996708, 'C', 'T', '', 0.04),
(6, 'chr1', 120069350, 120069350, 'G', 'A', '', 0.04),
(7, 'chr1', 120069420, 120069420, 'T', 'G', '', 0.04),
(8, 'chr1', 120069426, 120069426, '-', 'CCTCCTCCG', '', 0.14),
(9, 'chr5', 132589791, 132589791, 'G', 'C', '', 0.04),
(10, 'chr1', 13657, 13657, 'AG', '-', 'NA12878_38_Dragen_abra', 0.04),
(11, 'chr1', 17385, 17385, 'G', 'A', 'NA12878_38_Dragen_abra', 0.04),
(12, 'chr1', 494515, 494515, 'T', 'A', 'NA12878_38_Dragen_abra', 0.04),
(13, 'chr1', 826893, 826893, 'G', 'A', 'NA12878_38_Dragen_abra', 0.04),
(14, 'chr1', 827209, 827209, 'G', 'C', 'NA12878_38_Dragen_abra', 0.04),
(15, 'chr1', 827212, 827212, 'C', 'G', 'NA12878_38_Dragen_abra', 0.08),
(16, 'chr1', 827221, 827221, 'T', 'C', 'NA12878_38_Dragen_abra', 0.08),
(17, 'chr1', 827252, 827252, 'T', 'A', 'NA12878_38_Dragen_abra', 0.04),
(18, 'chr1', 931132, 931132, '-', 'CCCT', 'NA12878_38_Dragen_abra', 0.04),
(19, 'chr1', 935954, 935954, 'G', 'T', 'NA12878_38_Dragen_abra', 0.04),
(20, 'chr1', 941119, 941119, 'A', 'G', 'NA12878_38_Dragen_abra', 0.04),
(21, 'chr6', 81752209, 81752209, 'GACGCGCGGCGGCGGAGAGC', 'GCCCC', 'NA12878_38_Dragen_abra', 0.04);


--variant_classification
INSERT INTO variant_classification (id, variant_id, class, comment) VALUES
(1, 2, '5', 'pathogenic'),
(2, 8, '1', '');

--detected_variant
INSERT INTO detected_variant (processed_sample_id, variant_id, genotype) VALUES 
(1, 1, 'het'),
(1, 8, 'hom'),
(2, 1, 'het'),
(2, 3, 'het'),
(2, 4, 'het'),
(3, 1, 'hom'),
(3, 4, 'het'),
(4, 1, 'hom'),
(5, 10, 'het'),
(6, 10, 'het'),
(7, 10, 'het'),
(4, 11, 'het'),
(5, 11, 'het'),
(4, 12, 'het'),
(4, 21, 'het'),
(4, 14, 'hom'),
(5, 14, 'hom'),
(6, 14, 'hom'),
(4, 15, 'hom'),
(7, 14, 'hom'),
(5, 13, 'hom'),
(6, 21, 'hom'),
(7, 13, 'hom'),
(9, 13, 'hom'),
(11, 13, 'hom'),
(16, 13, 'hom'),
(10, 13, 'hom'),
(7, 4, 'hom'),
(10, 5, 'hom'),
(8, 10, 'hom'),
(4, 13, 'hom'),
(17, 13, 'hom');



--variant validation 
INSERT INTO variant_validation (id, user_id, sample_id, variant_id, genotype, status, comment) VALUES
(1, 2, 1, 8, 'hom', 'false positive', 'val com1'),
(2, 2, 2, 8, 'hom', 'false positive', 'val com1'),
(3, 2, 3, 8, 'hom', 'false positive', 'val com1'),
(4, 2, 4, 8, 'hom', 'false positive', 'val com1'),
(5, 2, 5, 8, 'hom', 'false positive', 'val com1'),
(6, 2, 4, 4, 'hom', 'true positive', 'val com2'),
(7, 2, 3, 4, 'hom', 'true positive', 'val com2'),
(8, 2, 2, 4, 'hom', 'true positive', 'val com2'),
(9, 2, 1, 4, 'hom', 'true positive', 'val com2'),
(10, 2, 5, 2, 'hom', 'true positive', 'val com2'),
(11, 2, 6, 7, 'hom', 'true positive', 'val com2'),
(12, 2, 2, 7, 'hom', 'true positive', 'val com2'),
(13, 2, 3, 7, 'hom', 'true positive', 'val com2'),
(14, 2, 5, 7, 'hom', 'true positive', 'val com2'),
(15, 2, 7, 7, 'hom', 'true positive', 'val com2'),
(16, 2, 4, 7, 'hom', 'true positive', 'val com2');

-- table `gene`
INSERT INTO `gene` (`id`, `hgnc_id`, `symbol`, `name`, `type`) VALUES
(1,1001, 'BRCA1','Breast cancer associated gene 1', 'protein-coding gene'),
(2,1002, 'BRCA2','Breast cancer associated gene 2', 'protein-coding gene'),
(3,1003, 'OR4F5', 'olfactory receptor family 4 subfamily F member 5', 'protein-coding gene'),
(4,1004, 'DIRC1', 'disrupted in renal carcinoma 1', 'protein-coding gene'),
(22712, 9121, 'PMS1', 'PMS1 homolog 1, mismatch repair system component', 'protein-coding gene');

-- table `gene_transcript`
INSERT INTO `gene_transcript` (`id`, `gene_id`, `name`, `source`, `chromosome`, `start_coding`, `end_coding`, `strand`) VALUES
(1, 1, 'uc001uua.1', 'ensembl', '13', 32899266, 32907523, '+'),
(2, 1, 'uc001uub.1', 'ensembl', '13', 32890598, 32972907, '+'),
(3, 1, 'uc031qky.1', 'ensembl', '13', 32929167, 32936796, '+'),
(4, 1, 'uc031qkz.1', 'ensembl', '13', null, null, '+'),
(5, 1, 'CCDS9344.1', 'ccds', '13', 32890598, 32973805, '+'),
(6, 2, 'ENST00000544455', 'ensembl', '13', 32889617, 32972907, '+'),
(7, 2, 'ENST00000380152', 'ensembl', '13', 32889611, 32973347, '+'),
(8, 3, 'ENST00000335137', 'ensembl', '1', 69091, 70008, '+'),
(9, 4, 'ENST00000308100', 'ensembl', '2', 189598882, 189654831, '-'),
(39236, 22712, 'uc010zfz.1', 'ensembl', '2', 190656536, 190670560, '+'),
(39237, 22712, 'uc010zga.1', 'ensembl', '2', 190656536, 190720597, '+'),
(39238, 22712, 'uc010zgb.1', 'ensembl', '2', 190656536, 190728952, '+'),
(39239, 22712, 'uc002urh.4', 'ensembl', '2', 190656536, 190742162, '+'),
(39240, 22712, 'uc002urk.4', 'ensembl', '2', 190656536, 190742162, '+'),
(39241, 22712, 'uc002uri.4', 'ensembl', '2', 190656536, 190742162, '+'),
(39242, 22712, 'uc010zgc.2', 'ensembl', '2', 190682853, 190742162, '+'),
(39243, 22712, 'uc010zgd.2', 'ensembl', '2', 190682853, 190742162, '+'),
(39244, 22712, 'uc002urj.3', 'ensembl', '2', null, null, '+'),
(39245, 22712, 'uc010fry.1', 'ensembl', '2', 190656536, 190738254, '+'),
(39246, 22712, 'uc010frz.3', 'ensembl', '2', 190656536, 190742162, '+'),
(39247, 22712, 'uc002url.3', 'ensembl', '2', 190682853, 190742162, '+'),
(39248, 22712, 'uc002urm.3', 'ensembl', '2', null, null, '+'),
(39249, 22712, 'uc002urn.1', 'ensembl', '2', 190718995, 190728841, '+'),
(85648, 22712, 'CCDS46474.1', 'ccds', '2', 190656536, 190742162, '+'),
(85649, 22712, 'CCDS46473.1', 'ccds', '2', 190656536, 190742162, '+'),
(85650, 22712, 'CCDS2302.1', 'ccds', '2', 190656536, 190742162, '+');

-- table `gene_exon`
INSERT INTO `gene_exon` (`transcript_id`, `start`, `end`) VALUES
(1, 32889617, 32889804),
(1, 32890559, 32890660),
(1, 32893214, 32893462),
(1, 32899213, 32899321),
(1, 32900238, 32900287),
(1, 32900379, 32900419),
(1, 32900636, 32900750),
(1, 32903580, 32903629),
(1, 32905056, 32905167),
(1, 32906409, 32907524),
(2, 32889617, 32889804),
(2, 32890559, 32890664),
(2, 32893214, 32893462),
(2, 32899213, 32899321),
(2, 32900238, 32900287),
(2, 32900379, 32900419),
(2, 32900636, 32900750),
(2, 32903580, 32903629),
(2, 32905056, 32905167),
(2, 32906409, 32907524),
(2, 32910402, 32915333),
(2, 32918695, 32918790),
(2, 32920964, 32921033),
(2, 32928998, 32929425),
(2, 32930565, 32930746),
(2, 32931879, 32932066),
(2, 32936660, 32936830),
(2, 32937316, 32937670),
(2, 32944539, 32944694),
(2, 32945093, 32945237),
(2, 32950807, 32950928),
(2, 32953454, 32953652),
(2, 32953887, 32954050),
(2, 32954144, 32954282),
(2, 32968826, 32969070),
(2, 32971035, 32971181),
(2, 32972299, 32973809),
(3, 32928998, 32929425),
(3, 32936660, 32936830),
(4, 32945093, 32945237),
(4, 32953454, 32953652),
(5, 32890598, 32890664),
(5, 32893214, 32893462),
(5, 32899213, 32899321),
(5, 32900238, 32900287),
(5, 32900379, 32900419),
(5, 32900636, 32900750),
(5, 32903580, 32903629),
(5, 32905056, 32905167),
(5, 32906409, 32907524),
(5, 32910402, 32915333),
(5, 32918695, 32918790),
(5, 32920964, 32921033),
(5, 32928998, 32929425),
(5, 32930565, 32930746),
(5, 32931879, 32932066),
(5, 32936660, 32936830),
(5, 32937316, 32937670),
(5, 32944539, 32944694),
(5, 32945093, 32945237),
(5, 32950807, 32950928),
(5, 32953454, 32953652),
(5, 32953887, 32954050),
(5, 32954144, 32954282),
(5, 32968826, 32969070),
(5, 32971035, 32971181),
(5, 32972299, 32972907),
(6, 32889617, 32889804),
(6, 32893214, 32893462),
(6, 32906409, 32907524),
(6, 32910402, 32915333),
(6, 32928998, 32929425),
(6, 32931879, 32932066),
(6, 32972299, 32973347),

(7, 32972299, 32973347),
(7, 32968826, 32969070),
(7, 32944539, 32944694),
(7, 32928998, 32929425),
(7, 32889611, 32889804),

(8, 69091, 70008),

(9, 189654590, 189654831),
(9, 189598882, 189599676),


(39236, 190648811, 190649319),
(39236, 190656516, 190656667),
(39236, 190660495, 190660677),
(39236, 190670378, 190670561),
(39237, 190648811, 190649319),
(39237, 190656516, 190656667),
(39237, 190660495, 190660677),
(39237, 190670378, 190670480),
(39237, 190682743, 190682906),
(39237, 190717381, 190717503),
(39237, 190718665, 190718808),
(39237, 190718965, 190719854),
(39237, 190720555, 190720642),
(39237, 190728469, 190728700),
(39238, 190648811, 190649319),
(39238, 190656516, 190656667),
(39238, 190670378, 190670480),
(39238, 190682743, 190682906),
(39238, 190708690, 190708806),
(39238, 190717381, 190717503),
(39238, 190718665, 190718808),
(39238, 190718965, 190719854),
(39238, 190728469, 190728954),
(39239, 190648811, 190649319),
(39239, 190656516, 190656667),
(39239, 190660495, 190660677),
(39239, 190670378, 190670480),
(39239, 190682743, 190682906),
(39239, 190708690, 190708806),
(39239, 190717381, 190717503),
(39239, 190718665, 190718808),
(39239, 190718965, 190719854),
(39239, 190728469, 190728954),
(39239, 190732525, 190732655),
(39239, 190738222, 190738382),
(39239, 190741998, 190742355),
(39240, 190648811, 190649319),
(39240, 190656516, 190656667),
(39240, 190660495, 190660677),
(39240, 190670378, 190670480),
(39240, 190682743, 190682906),
(39240, 190717381, 190717503),
(39240, 190718665, 190718808),
(39240, 190718965, 190719854),
(39240, 190728469, 190728954),
(39240, 190732525, 190732655),
(39240, 190738222, 190738382),
(39240, 190741998, 190742355),
(39241, 190648811, 190649319),
(39241, 190656516, 190656667),
(39241, 190660495, 190660677),
(39241, 190670378, 190670480),
(39241, 190682743, 190682906),
(39241, 190708690, 190708806),
(39241, 190717381, 190717503),
(39241, 190718665, 190718808),
(39241, 190718965, 190719854),
(39241, 190732525, 190732655),
(39241, 190738222, 190738382),
(39241, 190741998, 190742355),
(39242, 190648811, 190649319),
(39242, 190656516, 190656667),
(39242, 190660495, 190660677),
(39242, 190682743, 190682906),
(39242, 190708690, 190708806),
(39242, 190717381, 190717503),
(39242, 190718665, 190718808),
(39242, 190718965, 190719854),
(39242, 190728469, 190728954),
(39242, 190732525, 190732655),
(39242, 190738222, 190738382),
(39242, 190741998, 190742355),
(39243, 190648811, 190649319),
(39243, 190656516, 190656667),
(39243, 190682743, 190682906),
(39243, 190708690, 190708806),
(39243, 190717381, 190717503),
(39243, 190718665, 190718808),
(39243, 190718965, 190719854),
(39243, 190728469, 190728954),
(39243, 190732525, 190732655),
(39243, 190738222, 190738382),
(39243, 190741998, 190742355),
(39244, 190649239, 190649290),
(39244, 190650072, 190650197),
(39244, 190656516, 190656667),
(39244, 190660495, 190660677),
(39244, 190670378, 190670480),
(39244, 190682743, 190682906),
(39244, 190717381, 190717503),
(39244, 190718665, 190718808),
(39244, 190718965, 190719854),
(39244, 190728469, 190728700),
(39244, 190732525, 190732655),
(39244, 190738222, 190738382),
(39244, 190741998, 190742355),
(39245, 190656516, 190656667),
(39245, 190660495, 190660677),
(39245, 190670378, 190670480),
(39245, 190682743, 190682906),
(39245, 190717381, 190717503),
(39245, 190718665, 190718808),
(39245, 190718965, 190719854),
(39245, 190728469, 190728700),
(39245, 190738222, 190738382),
(39246, 190656516, 190656667),
(39246, 190660495, 190660677),
(39246, 190670378, 190670480),
(39246, 190682743, 190682906),
(39246, 190741998, 190742355),
(39247, 190682743, 190682906),
(39247, 190717381, 190717503),
(39247, 190718665, 190718808),
(39247, 190718965, 190719854),
(39247, 190732525, 190732655),
(39247, 190738222, 190738382),
(39247, 190741998, 190742355),
(39248, 190708690, 190708806),
(39248, 190717381, 190717503),
(39248, 190718665, 190718808),
(39248, 190718965, 190719854),
(39248, 190728469, 190728841),
(39248, 190732525, 190732655),
(39248, 190738222, 190738382),
(39248, 190741998, 190742355),
(39249, 190718655, 190718808),
(39249, 190718965, 190719854),
(39249, 190728469, 190728841),
(85648, 190656536, 190656667),
(85648, 190660495, 190660677),
(85648, 190670378, 190670480),
(85648, 190682743, 190682906),
(85648, 190717381, 190717503),
(85648, 190718665, 190718808),
(85648, 190718965, 190719854),
(85648, 190728469, 190728954),
(85648, 190732525, 190732655),
(85648, 190738222, 190738382),
(85648, 190741998, 190742162),
(85649, 190656536, 190656667),
(85649, 190660495, 190660677),
(85649, 190670378, 190670480),
(85649, 190682743, 190682906),
(85649, 190708690, 190708806),
(85649, 190717381, 190717503),
(85649, 190718665, 190718808),
(85649, 190718965, 190719854),
(85649, 190732525, 190732655),
(85649, 190738222, 190738382),
(85649, 190741998, 190742162),
(85650, 190656536, 190656667),
(85650, 190660495, 190660677),
(85650, 190670378, 190670480),
(85650, 190682743, 190682906),
(85650, 190708690, 190708806),
(85650, 190717381, 190717503),
(85650, 190718665, 190718808),
(85650, 190718965, 190719854),
(85650, 190728469, 190728954),
(85650, 190732525, 190732655),
(85650, 190738222, 190738382),
(85650, 190741998, 190742162);

INSERT INTO `geneinfo_germline`(`symbol`, `inheritance`, `gnomad_oe_mis`, `gnomad_oe_syn`, `gnomad_oe_lof`, `comments`) VALUES
('BRCA1', 'AD', 0.00, 0.00, 0.00, ''),
('BRCA2', 'AD', NULL, NULL, NULL, '');

INSERT INTO `hpo_term` (`id`, `hpo_id`, `name`, `definition`, `synonyms`) VALUES
(1, 'HP:0000007', 'Autosomal recessive inheritance', '\"A mode of inheritance that is observed for traits related to a gene encoded on one of the autosomes (i.e., the human chromosomes 1-22) in which a trait manifests in homozygotes. In the context of medical genetics, autosomal recessive disorders manifest in homozygotes (with two copies of the mutant allele) or compound heterozygotes (whereby each copy of a gene has a distinct mutant allele).', ''),
(2, 'HP:0001427', 'Mitochondrial inheritance', '\"A mode of inheritance that is observed for traits related to a gene encoded on the mitochondrial genome. Because the mitochondrial genome is essentially always maternally inherited, a mitochondrial condition can only be transmitted by females, although the condition can affect both sexes. The proportion of mutant mitochondria can vary (heteroplasmy).', '');

INSERT INTO `hpo_genes` (`hpo_term_id`, `gene`, `details`, `evidence`) VALUES
(1, 'OR4F5', '(ClinVar,,n/a); (HPO,,n/a)', 'n/a'),
(2, 'OR4F5', '(OMIM,(2), low)', 'low');

