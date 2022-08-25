import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
import json
import re
import argparse
import random
import subprocess
import multiprocessing as mp
from multiprocessing import freeze_support
from os.path import exists

DESCRIPTION = "Runs blastn for a fasta of predicted viral contigs agains reference viruses and HERVs."
contaminators = ['AC_000018','AC_000007','AC_000008','AC_000017','NC_001405','NC_001460','NC_001454','NC_010956','NC_011202','NC_011203','NC_012959','AF192278','AB100179','AY884927','AB100229','AF017224','JX012349','AB100191','AB811831','AB120904','KX814024','U07813','KX814024','AB022748','AB365180','AY523480','U88895','M87811','AB100179','AB100179','AK037062','AF182203','AF229251','AJ290702','D12848','EF102092','KY865072','KY865082','U48716','GU939858','DQ358089','M21123','M80675','X17124','X51459','X94355','AF017203','AF050455', 'AF464738', 'AH015066', 'AF059603','AF059603','EU727812','GU940037','GU940175','AJ304824', 'AM229312', 'AY533034', 'AY864614', 'EF133960', 'EF468506', 'FJ600523', 'GQ326556', 'HM131064', 'JN415766', 'KX759202', 'KY910034','GU939840','D90663', 'DQ536167', 'DQ985805', 'EF486298', 'EU730688', 'JQ726570', 'Y18778', 'JQ790991','AB022148','Y18782', 'Z80017', 'Z80031','AF378027', 'LN879542','AF378029', 'AY204209', 'AY641409', 'D90698', 'EF469548', 'JX997156','NC_009503','Y18778','KF214645','AF182203','EF102092','KY865072','KY865077','AJ965724','AF049110', 'AF050439', 'AF090446', 'AF464738', 'AF464773', 'AY234827', 'EF468508', 'EF562447', 'KT724926', 'KX759202', 'MH675928', 'U41000', 'U68408','EU727822','GU940037','GU940473','KY684104','NC_009503','KU978791','NC_026440','AF216200','KU978791','AF216200','AJ290690','EF102092','JQ975222','KY865039','KY865055','KY865072','AB021264','KX857216', 'AB100177', 'AB120923', 'HQ536016','KF729698','KF729762','CY125728','FJ610152', 'AB079295','M87813', 'AB120932', 'KF729762', 'AB071970', 'AF118055', 'KC616429','AY299398','AF315089','V00040','NC_040309','AB008772', 'AF396663', 'AF453669', 'AJ270051', 'AJ582609', 'AY641408', 'D90603', 'D90632', 'FN432151', 'JN135383', 'JN800031', 'DQ011795', 'GU230478', 'DQ011795', 'AF104027', 'KF198064', 'DQ011795', 'DQ011795','KR873417', 'KY404253', 'KY404254', 'Y18777', 'Z80002', 'Z84554', 'AJ291362', 'AB023927', 'AB022149', 'EU309022','NC_005831','MH155897','NC_020810','NC_038603','NC_044855','KU133671','KF892052','AB022149','KC616429','MH678757','AF315098','AB811832','D14335','AH001391','HQ536016','AB079295','AB071970','AF017189','AF017187','D14335','NC_044046','AF001600','AF050215','AF065434','AF069450','AJ251464','FM212573','KY684103','KY910011','AF106947','GU934326','AJ291362','NC_002645','NC_035889','AH001391','NC_038496','D55717','FJ197996','U96004','MG451071','AB055953', 'AB022147', 'AB022150', 'AB022153', 'AB022154', 'AB022745', 'AB022777', 'AB023512', 'AB023514', 'AB023517', 'AB023535', 'AB023541', 'AB023687', 'AB023688', 'AB023690', 'AB023691', 'AB023692', 'AB023694', 'AB023695', 'AB023696', 'AB023697', 'AB023698', 'AB023699', 'AB023700', 'AB023702', 'AB023703', 'AB023704', 'AB023705', 'AB023706', 'AB023707', 'AB023708', 'AB023710', 'AB023711', 'AB023914', 'AB023916', 'AB023917', 'AB023931', 'AB023933', 'AB025947', 'AB025948', 'AB025949', 'AB025950', 'AB025951', 'AB025952', 'AB025953', 'AB025954', 'AB025955', 'AB025957', 'AB025960', 'AB025963', 'AB047209', 'AB047240', 'AB048243', 'AB052104', 'AB055951', 'AB059366', 'AB059372', 'AB060008', 'AB060026', 'AB060027', 'AB060039', 'AB062377', 'AB062378', 'AB062379', 'AB100197', 'AB100207', 'AB100212', 'AB100219', 'AB100220', 'AB100233', 'AB111966', 'AB111967', 'AB120307', 'AB120770', 'AB120895', 'AB120902', 'AB120922', 'AB164262', 'AB164264', 'AB168289', 'AB185332', 'AB214978', 'AB439574', 'AB443932', 'AB447373', 'AB447374', 'AB453920', 'AB728590', 'AB749815', 'AB811796', 'AB811797', 'AB811801', 'AB811803', 'AB811804', 'AB811806', 'AB811808', 'AB811810', 'AB811813', 'AB811814', 'AB811821', 'AB811825', 'AB811827', 'AB811828', 'AB811834', 'AB811835', 'AB811836', 'AB811837', 'AB811838', 'AB811839', 'AB811840', 'AB811841', 'AB811842', 'AB811844', 'AB897842', 'AB972431', 'AC146851', 'AF012331', 'AF012335', 'AF015685', 'AF017192', 'AF017193', 'AF017209', 'AF017213', 'AF017214', 'AF017215', 'AF017218', 'AF017221', 'AF017222', 'AF017225', 'AF017226', 'AF017227', 'AF026248', 'AF026249', 'AF026250', 'AF026251', 'AF026252', 'AF026253', 'AF026254', 'AF026255', 'AF064191', 'AF065659', 'AF065660', 'AF065684', 'AF065688', 'AF065689', 'AF065698', 'AF065699', 'AF065702', 'AF065724', 'AF065725', 'AF065755', 'AF065756', 'AF067487', 'AF074086', 'AF078838', 'AF083247', 'AF087913', 'AF094515', 'AF104019', 'AF104020', 'AF104021', 'AF104022', 'AF104023', 'AF104024', 'AF104025', 'AF104026', 'AF104028', 'AF104029', 'AF110315', 'AF126162', 'AF134163', 'AF134164', 'AF135486', 'AF139840', 'AF141972', 'AF164611', 'AF164615', 'AF172456', 'AF172459', 'AF172490', 'AF172491', 'AF172492', 'AF172494', 'AF172496', 'AF172497', 'AF172498', 'AF172500', 'AF180430', 'AF182201', 'AF191073', 'AF193443', 'AF216197', 'AF216198', 'AF216201', 'AF216202', 'AF227993', 'AF227994', 'AF227995', 'AF227997', 'AF227998', 'AF228933', 'AF229252', 'AF244571', 'AF260248', 'AF260251', 'AF260253', 'AF268174', 'AF275345', 'AF288220', 'AF288641', 'AF290423', 'AF294261', 'AF294264', 'AF311937', 'AF315093', 'AF315094', 'AF315100', 'AF330705', 'AF333071', 'AF333072', 'AF339051', 'AF378030', 'AF378031', 'AF378032', 'AF378034', 'AF378035', 'AF380292', 'AF387849', 'AF390031', 'AF394944', 'AF396664', 'AF397035', 'AF411058', 'AF411059', 'AF440570', 'AF450249', 'AF453664', 'AF453665', 'AF453674', 'AF453676', 'AF453677', 'AF453679', 'AF453680', 'AF453683', 'AF453684', 'AF499232', 'AF531472', 'AH001392', 'AH001458', 'AH001459', 'AH001476', 'AH002364', 'AH005269', 'AH007031', 'AH007766', 'AH007767', 'AH008411', 'AH008413', 'AH012932', 'AH012974', 'AH013886', 'AH015602', 'AJ000887', 'AJ000888', 'AJ000889', 'AJ000890', 'AJ000891', 'AJ000892', 'AJ001235', 'AJ002617', 'AJ002619', 'AJ004945', 'AJ004947', 'AJ005412', 'AJ131152', 'AJ131153', 'AJ131154', 'AJ131155', 'AJ131156', 'AJ131157', 'AJ131158', 'AJ131159', 'AJ131160', 'AJ223973', 'AJ244585', 'AJ269530', 'AJ279072', 'AJ288180', 'AJ289709', 'AJ289710', 'AJ289711', 'AJ290613', 'AJ291349', 'AJ291351', 'AJ291352', 'AJ291353', 'AJ291356', 'AJ291357', 'AJ291359', 'AJ291360', 'AJ291716', 'AJ291759', 'AJ293237', 'AJ301642', 'AJ303092', 'AJ304461', 'AJ421984', 'AJ421985', 'AJ421986', 'AJ421987', 'AJ421989', 'AJ421990', 'AJ575660', 'AJ582608', 'AJ582610', 'AJ582611', 'AJ582612', 'AJ879490', 'AK022175', 'AK024677', 'AK027132', 'AK027691', 'AK055945', 'AK074464', 'AK091794', 'AK092197', 'AK093935', 'AK122924', 'AM698090', 'AM901569', 'AM909712', 'AM909713', 'AM909714', 'AP018495', 'AY015682', 'AY029538', 'AY032741', 'AY034063', 'AY034066', 'AY040832', 'AY040833', 'AY054379', 'AY178976', 'AY204210', 'AY224336', 'AY224342', 'AY224368', 'AY224376', 'AY224377', 'AY299397', 'AY371031', 'AY395526', 'AY397620', 'AY423624', 'AY425963', 'AY513268', 'AY525758', 'AY535425', 'AY573559', 'AY573561', 'AY641410', 'AY650289', 'AY691949', 'AY747601', 'AY846360', 'AY878356', 'AY884838', 'AY884840', 'AY884841', 'AY884844', 'AY884845', 'AY884847', 'AY884859', 'AY884865', 'AY884891', 'AY884893', 'AY884894', 'AY884895', 'AY884899', 'AY884907', 'AY884911', 'AY884913', 'AY884928', 'AY884931', 'AY884934', 'AY884936', 'AY884939', 'AY925147', 'AY925148', 'AY945804', 'AY945805', 'BC011670', 'BC012129', 'BC018736', 'BC032544', 'BC033877', 'BC151135', 'BK000394', 'BK009405', 'D00735', 'D00835', 'D10083', 'D10450', 'D11078', 'D12832', 'D87055', 'D87056', 'D87057', 'D87058', 'D90606', 'D90608', 'D90629', 'D90639', 'D90646', 'D90648', 'D90671', 'D90677', 'D90680', 'D90682', 'D90683', 'D90696', 'DQ112093', 'DQ112094', 'DQ112095', 'DQ112096', 'DQ112097', 'DQ112098', 'DQ112101', 'DQ112104', 'DQ112105', 'DQ112108', 'DQ112111', 'DQ112119', 'DQ112122', 'DQ112123', 'DQ112126', 'DQ112127', 'DQ112128', 'DQ112130', 'DQ112131', 'DQ112132', 'DQ112134', 'DQ112135', 'DQ112136', 'DQ112137', 'DQ112141', 'DQ112142', 'DQ112143', 'DQ112144', 'DQ112148', 'DQ112149', 'DQ112150', 'DQ112152', 'DQ112154', 'DQ118011', 'DQ127584', 'DQ127741', 'DQ127777', 'DQ127791', 'DQ128036', 'DQ157732', 'DQ273234', 'DQ273251', 'DQ273259', 'DQ273262', 'DQ273264', 'DQ360593', 'DQ445232', 'DQ445236', 'DQ445241', 'DQ445243', 'DQ445248', 'DQ445253', 'DQ445619', 'DQ445623', 'DQ445624', 'DQ473744', 'DQ666286', 'DQ838494', 'DQ857990', 'DQ983234', 'DQ985749', 'DQ985751', 'DQ985755', 'DQ985766', 'DQ985767', 'DQ985801', 'EF066511', 'EF102090', 'EF179138', 'EF194101', 'EF453756', 'EF576506', 'EF629313', 'EU001284', 'EU009623', 'EU009624', 'EU105454', 'EU105459', 'EU195798', 'EU306659', 'EU308984', 'EU309011', 'EU410304', 'EU646311', 'EU646322', 'EU669866', 'EU727925', 'EU727933', 'EU728030', 'EU728134', 'EU791617', 'EU812120', 'EU862277', 'EU963184', 'FJ173819', 'FJ173823', 'FJ435168', 'FJ435169', 'FJ539167', 'FJ593289', 'FJ640151', 'FJ751933', 'FJ788422', 'FJ796235', 'FM210528', 'FM210529', 'FM212574', 'FN432141', 'FN432142', 'FN432143', 'FN432144', 'FN432145', 'FN432149', 'FN432150', 'GQ183504', 'GQ252855', 'GQ252887', 'GQ267466', 'GQ294555', 'GU317962', 'GU476554', 'GU592189', 'GU939763', 'GU940376', 'GU940622', 'GU940623', 'GU980198', 'HE984505', 'HE984515', 'HE984519', 'HE984522', 'HE984542', 'HE984564', 'HE984565', 'HF548798', 'HF559481', 'HG523434', 'HG531765', 'HM124452', 'HM133903', 'HM454274', 'HM747958', 'HM849633', 'HM991350', 'HQ442266', 'HQ457130', 'HQ540592', 'HQ540594', 'HQ910537', 'JF430693', 'JF430694', 'JF513052', 'JN185000', 'JN185001', 'JN185002', 'JN185003', 'JN204354', 'JN379815', 'JN624480', 'JN675007', 'JN675008', 'JN675010', 'JN675012', 'JN675016', 'JN675018', 'JN675019', 'JN675020', 'JN675022', 'JN675028', 'JN675030', 'JN675031', 'JN675032', 'JN675038', 'JN675045', 'JN675047', 'JN675053', 'JN675054', 'JN675055', 'JN675056', 'JN675057', 'JN675059', 'JN675060', 'JN675062', 'JN675065', 'JN675070', 'JN675071', 'JN675078', 'JN675079', 'JN675086', 'JN675088', 'JN675094', 'JN817430', 'JN969602', 'JN969603', 'JQ043204', 'JQ314449', 'JQ354986', 'JQ726569', 'JQ790928', 'JQ790976', 'JQ790982', 'JQ790990', 'JQ811903', 'JQ966584', 'JQ966586', 'JQ966587', 'JQ966590', 'JQ966591', 'JQ975179', 'JQ975186', 'JQ975233', 'JQ999963', 'JX120155', 'JX120157', 'JX120159', 'JX185662', 'JX185663', 'JX185666', 'JX185667', 'JX291540', 'JX904130', 'JX904314', 'JX997155', 'JX997185', 'JX997187', 'K01001', 'K01100', 'K01103', 'K02533', 'KC440852', 'KC503852', 'KC519456', 'KC519458', 'KC611152', 'KC786228', 'KF254348', 'KF254380', 'KF254387', 'KF254392', 'KF261120', 'KF493876', 'KF493877', 'KF651980', 'KF779381', 'KF892040', 'KF898354', 'KF913885', 'KJ000234', 'KJ191556', 'KJ410306', 'KJ716849', 'KJ946251', 'KM043884', 'KM254164', 'KM254174', 'KM577554', 'KP162058', 'KP183324', 'KP183325', 'KP183327', 'KP258515', 'KP317931', 'KR029591', 'KR029607', 'KR873406', 'KR873412', 'KR873415', 'KT180230', 'KT203716', 'KT427527', 'KT427528', 'KT715968', 'KT724897', 'KT724934', 'KT895861', 'KT895870', 'KT982173', 'KT982185', 'KU052576', 'KU054256', 'KU054260', 'KU054266', 'KU054268', 'KU054271', 'KU576450', 'KU576959', 'KU746283', 'KU758926', 'KX018009', 'KX689266', 'KX759199', 'KX759200', 'KX775983', 'KX775984', 'KX775986', 'KX775988', 'KX775990', 'KX775991', 'KX775992', 'KX775993', 'KX775994', 'KX775995', 'KX775996', 'KX775997', 'KX775998', 'KX813893', 'KX814186', 'KX814288', 'KX814289', 'KX814290', 'KX814299', 'KX814301', 'KX814302', 'KX814303', 'KX814304', 'KX814313', 'KX814314', 'KX853510', 'KY052815', 'KY052857', 'KY404242', 'KY404249', 'KY404250', 'KY404252', 'KY404265', 'KY404271', 'KY523104', 'KY559403', 'KY684114', 'KY684123', 'KY766069', 'KY852349', 'KY865041', 'KY865051', 'KY865078', 'KY865083', 'KY865088', 'KY883317', 'KY883335', 'L13783', 'L19439', 'L22023', 'L40326', 'LC012636', 'LC012641', 'LC176792', 'LC176796', 'LC193725', 'LC315647', 'LN680393', 'LR584255', 'LS992247', 'LT837585', 'LT907991', 'M10976', 'M18048', 'M18706', 'M19755', 'M20211', 'M23367', 'M27088', 'M27826', 'M27828', 'M30839', 'M33842', 'M34549', 'M57950', 'M59893', 'M59901', 'M62431', 'M73484', 'M85292', 'M87812', 'M92449', 'MF069043', 'MF094128', 'MF094130', 'MF405918', 'MF468142', 'MF468144', 'MF468147', 'MF624891', 'MF624974', 'MF625012', 'MF770979', 'MG298823', 'MG298860', 'MG298870', 'MG598800', 'MG598801', 'MG598802', 'MG679361', 'MG770349', 'MG770354', 'MG934417', 'MG986233', 'MG986280', 'MH052023', 'MH057869', 'MH057876', 'MH057879', 'MH057888', 'MH057893', 'MH057902', 'MH057907', 'MH057919', 'MH057925', 'MH057934', 'MH057950', 'MH057951', 'MH057954', 'MH057959', 'MH057962', 'MH057965', 'MH057990', 'MH057992', 'MH057993', 'MH078547', 'MH121119', 'MH319740', 'MH319764', 'MH552529', 'MH558113', 'MH590376', 'MH590398', 'MH590409', 'MH590442', 'MH590571', 'MH645153', 'MH645156', 'MH678754', 'MH678755', 'MH678759', 'MH678766', 'MH678768', 'MH678771', 'MH678773', 'MH678778', 'MH678779', 'MH678781', 'MH678785', 'MH678787', 'MH771723', 'MH883318', 'MH892403', 'MK039127', 'MK071981', 'MK100570', 'MK458165', 'MK500307', 'MK500336', 'MN369532', 'MN646695', 'MN887110', 'MT222957', 'NC_001422', 'NC_001798', 'NC_008168', 'NC_012986', 'NC_013756', 'NC_014372', 'NC_018464', 'NC_018476', 'NC_022518', 'NC_024450', 'NC_024625', 'NC_028045', 'NC_030230', 'NC_032001', 'NC_032111', 'NC_033774', 'NC_035758', 'NC_038727', 'NC_040306', 'NC_055230', 'NC_055235', 'S46400', 'S46401', 'S61070', 'U07814', 'U09236', 'U15647', 'U23467', 'U25276', 'U27612', 'U27973', 'U29659', 'U34991', 'U37066', 'U47118', 'U76261', 'U86599', 'U86698', 'U88896', 'U88901', 'U92818', 'U92819', 'U92820', 'U95997', 'U95998', 'U95999', 'U96000', 'U96001', 'U96002', 'U96003', 'U96006', 'U96011', 'U96043', 'U96044', 'U96046', 'U96047', 'U96049', 'U96050', 'U96052', 'U96056', 'U96058', 'U96061', 'U96065', 'U96070', 'U96072', 'V01571', 'V01573', 'X00255', 'X01172', 'X05349', 'X06275', 'X12716', 'X12719', 'X12720', 'X12721', 'X12722', 'X12723', 'X12724', 'X14949', 'X14950', 'X56464', 'X57147', 'X57168', 'X59357', 'X63184', 'X67284', 'X72790', 'X82271', 'X83497', 'X87297', 'X94741', 'X96883', 'Y17833', 'Y18773', 'Y18779', 'Y18781', 'Y18783', 'Z14310', 'Z17327', 'Z21850', 'Z21851', 'Z21852', 'Z25525', 'Z48163', 'Z48633', 'Z49885', 'Z54175', 'Z72499', 'Z72519', 'Z80000', 'Z80007', 'Z80009', 'Z80022', 'Z80026', 'Z80027', 'Z80030', 'Z80033', 'Z80035', 'Z80038', 'Z80040', 'Z80043', 'Z80044', 'Z80045', 'Z80047', 'Z80049', 'Z80064', 'Z80072', 'Z80075', 'Z84550', 'Z95336', 'Z95340']

words_rm = ['sapiens','ORF','clone','Vaccinia','endogenous','UNVERIFIED','gene','Shamonda','Tokyovirus','BeAn','Pepper',
            'mRNA','Human DNA','retrotransposon','partial','HIV-1 isolate','UTR','viroid','ERV','RTVL','Towne','White spot','Uncultured',
            'long terminal repeat','cell line','ntegration site','Equid','Suid','Semliki forest','defective','segment','CRESS','cRNA',
            'Choristoneura','DNA','noncoding','Endogenous','Saccharomyces','Macacine betaherpesvirus','RNA','Macacine betaherpesvirus',
            'IgG','Acanthocystis','Klosneuvirus','Pithovirus','mosaic','transposon','genome assembly','retroelement','LTR','Ullucus tymovirus']
words_rm=[]

TCGA=False
#https://hive.biochemistry.gwu.edu/rvdb/
gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

basepairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def handle_non_ATGC(sequence):
    """
    Handle non ATGCs.
    :param sequence: String input.
    :return: String output (only ATCGs), with randomly assigned bp to non-ATGCs.
    """
    DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
    NUCLEOTIDES = sorted([x for x in DEFAULT_NUC_ORDER.keys()])

    ret = re.sub('[^ATCG]', random.choice(NUCLEOTIDES), sequence)
    assert len(ret) == len(sequence)
    return ret

def translate_frameshifted(sequence, gcode=gencode):
    """
    Translation of nucleotide to amino acid
    :param sequence: a section of nucleotide sequence
    :param gcode: gencode dictionary
    :return: amino acide sequence
    """
    translate = ''.join([gcode.get(sequence[3 * i:3 * i + 3]) for i in range(len(sequence) // 3)])
    return translate

def reverse_complement(sequence, bpairs=basepairs):
    """
    Genertate the reversed sequence
    :param sequence: a section of nucleotide sequence
    :param bpairs: basepairs dictionary
    :return: the reversed string of the nucleotide sequence
    """
    reversed_sequence = (sequence[::-1])
    rc = ''.join([bpairs.get(reversed_sequence[i]) for i in range(len(sequence))])
    return rc

def get_project_id(smp_name):

    if TCGA:
        df_sample = pd.read_csv('gdc_sample.csv')

        filename = list(df_sample['File Name'])
        f = [i.replace('_gdc_realn_rehead.bam', '') for i in filename]
        caseid = list(df_sample['Case ID'])
        sampletp = list(df_sample['Sample Type'])
        projectid = list(df_sample['Project ID'])
        findi = lambda x, v: [i for i, j in enumerate(x) if j == v]
        case_id = []
        sample_tp = []
        project_id = []
        for i in range(len(smp_name)):

            indexes = [k for k in findi(f, smp_name[i])]
            c_id = [caseid[k] for k in indexes][0]
            s_tp = [sampletp[k] for k in indexes][0]
            p_id = [projectid[k] for k in indexes][0]
            case_id.append(c_id)
            sample_tp.append(s_tp)
            project_id.append(p_id)
    else:
        case_id = smp_name
        sample_tp = smp_name
        project_id = smp_name

    return case_id, sample_tp, project_id

def get_table_reference(infile, outfie, filen):
    '''

    :param infile: input file (contigs-fasta), at dbgap-mutporg/
    :param outfie: output file, at blast
    :param filen: generated files name (dirname_ref+'/'+file.split('/')[-1].split('.')[0])
    :return: combined csv data with refseq families

    '''
    findi = lambda x, v: [i for i, j in enumerate(x) if j == v]
    get_value = lambda x, y: [x.get(i) for i in y]

    col_names = ["id", "virus", "score", "coverage", "b", "c", "d", "e", "f", "g", "h", "i"]
    out11 = pd.read_csv(outfie, sep='\t', header=None, names=col_names)

    fasta_ids = list(out11['id'])
    fasta_id = [i.split('_bact')[0] for i in fasta_ids]

    virs = list(out11['virus'])
    vir = [i.split('|')[0].split('.')[0] for i in virs]

    uni_vir = sorted(list(set(vir)))
    score = list(out11['score'])
    coverage = list(out11['coverage'])

    ref_dict = {'NC_055230': 'Akhmeta virus isolate Akhmeta_2013-88, complete genome', 'NC_055330': 'Alenquer virus segment L, complete sequence', 'NC_055331': 'Alenquer virus segment S, complete sequence', 'NC_055332': 'Alenquer virus segment M, complete sequence', 'NC_055339': 'Echarte virus segment L, complete sequence', 'NC_055340': 'Echarte virus segment S, complete sequence', 'NC_055341': 'Echarte virus segment M, complete sequence', 'NC_055342': 'Maldonado virus strain FMD 0077 segment S nonstructural protein and nucleocapsid genes, complete cds', 'NC_055343': 'Maldonado virus strain FMD 0077 segment M polyprotein gene, complete cds', 'NC_055344': 'Maldonado virus strain FMD 0077 segment L L protein gene, complete cds', 'NC_045512': 'Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome', 'NC_044932': 'Norovirus GII GII.NA2[PNA2], complete sequence', 'NC_044853': 'Norovirus GI strain Hu/JP/1998/GI.6[PNA4]/No20-Saitama-98-17, complete genome', 'NC_044854': 'Norovirus GI strain Hu/JP/2000/GI.6[PNA1]/WUG1, complete genome', 'NC_044855': 'Norovirus GIV strain Hu/US/2016/GIV.NA1[PNA1]/WI7002, complete genome', 'NC_044856': 'Norovirus GI strain Hu/BD/2011/GI.7[PNA2]/Dhaka1882, complete genome', 'NC_044045': 'Norovirus GII/Hu/JP/2007/GII.P15_GII.15/Sapporo/HK299, complete genome', 'NC_044046': 'Norovirus GII/Hu/JP/2011/GII/Yuzawa/Gira2HS, complete genome', 'NC_043445': 'Guenon simian foamy virus isolate AG16, complete genome', 'NC_043450': 'Severe fever with thrombocytopenia syndrome virus strain HNXH segment L, complete sequence', 'NC_043451': 'Severe fever with thrombocytopenia syndrome virus strain HNXH segment M, complete sequence', 'NC_043452': 'Severe fever with thrombocytopenia syndrome virus strain HNXH segment S, complete sequence', 'NC_043585': 'Ilesha virus strain R5964 segment S, complete sequence', 'NC_043586': 'Ilesha virus strain R5964 segment M, complete sequence', 'NC_043587': 'Ilesha virus strain R5964 segment L, complete sequence', 'NC_043615': 'Fort Sherman virus strain 86MSP18 segment M, complete sequence', 'NC_043616': 'Fort Sherman virus strain 86MSP18 segment S, complete sequence', 'NC_043617': 'Fort Sherman virus strain 86MSP18 segment L, complete sequence', 'NC_040876': 'Norovirus GII.P7_GII.6, complete genome', 'NC_040550': 'Gammapapillomavirus sp. isolate Gamma_w34c04a, complete genome', 'NC_040619': 'Gammapapillomavirus sp. isolate Gamma_LCOSOc196, complete genome', 'NC_040620': 'Gammapapillomavirus sp. isolate GammaDyskD_w07c34d, complete genome', 'NC_040688': 'Gammapapillomavirus sp. isolate Gamma_w11C51, complete genome', 'NC_040691': 'Gammapapillomavirus sp. isolate Gammai915_ga2c70, complete genome', 'NC_040803': 'Gammapapillomavirus sp. isolate Gamma_w23c08c, complete genome', 'NC_040804': 'Gammapapillomavirus sp. isolate Gamma_w27c157c, complete genome', 'NC_040805': 'Gammapapillomavirus sp. isolate GammaiCH2_EV03c434, complete genome', 'NC_040806': 'Gammapapillomavirus sp. isolate Gamma_w27c03a, complete genome', 'NC_039897': 'Norovirus GI/Hu/JP/2007/GI.P3_GI.3/Shimizu/KK2866, complete genome', 'NC_039475': 'Norovirus GII.17, complete genome', 'NC_039477': 'Norovirus GII, complete genome', 'NC_038283': 'Ekpoma virus 2 isolate EKV-2, partial genome', 'NC_038294': 'Betacoronavirus England 1 isolate H123990006, complete genome', 'NC_038307': 'Coxsackievirus B3 mRNA, complete genome', 'NC_038337': 'Torque teno virus 9 clone BM1C-18 ORF2 (ORF2) gene, complete cds; and nonfunctional ORF1 (ORF1) gene, complete sequence', 'NC_038345': 'Torque teno mini virus 10 isolate LIL-y1 ORF2, ORF1, ORF3, and ORF4 genes, complete cds', 'NC_038346': 'Torque teno mini virus 11 isolate LIL-y2 ORF2, ORF1, and ORF3 genes, complete cds', 'NC_038347': 'Torque teno mini virus 12 isolate LIL-y3 ORF2, ORF1, ORF3, and ORF4 genes, complete cds', 'NC_038350': 'Torque teno midi virus 3 isolate 2PoSMA ORF2 and ORF1 genes, complete cds', 'NC_038351': 'Torque teno midi virus 4 isolate 6PoSMA ORF2, ORF1, and ORF3 genes, complete cds', 'NC_038352': 'Torque teno midi virus 5 DNA, complete genome, isolate: MDJHem2', 'NC_038353': 'Torque teno midi virus 6 DNA, complete genome, isolate: MDJHem3-1', 'NC_038354': 'Torque teno midi virus 7 DNA, complete genome, isolate: MDJHem3-2', 'NC_038355': 'Torque teno midi virus 8 DNA, complete genome, isolate: MDJN1', 'NC_038356': 'Torque teno midi virus 9 DNA, complete genome, isolate: MDJN2', 'NC_038357': 'Torque teno midi virus 10 DNA, complete genome, isolate: MDJN14', 'NC_038358': 'Torque teno midi virus 11 DNA, complete genome, isolate: MDJN47', 'NC_038359': 'Torque teno midi virus 12 DNA, complete genome, isolate: MDJN51', 'NC_038360': 'Torque teno midi virus 13 DNA, complete genome, isolate: MDJN69', 'NC_038361': 'Torque teno midi virus 14 DNA, complete genome, isolate: MDJN97', 'NC_038412': 'Cyclovirus PK5510, complete genome', 'NC_038413': 'Cyclovirus PK5222, complete genome', 'NC_038414': 'Cyclovirus TN25, complete genome', 'NC_038415': 'Cyclovirus PK5034, complete genome', 'NC_038416': 'Cyclovirus NG12, complete genome', 'NC_038417': 'Cyclovirus NG14, complete genome', 'NC_038418': 'Cyclovirus SL-108277, complete genome', 'NC_038426': 'Hepacivirus B polypeptide gene, complete cds', 'NC_038496': 'Gemycircularvirus C1c capsid protein (Cap) gene, complete cds', 'NC_038594': 'Lebombo virus isolate SAAR 3896 segment 1, complete sequence', 'NC_038595': 'Lebombo virus isolate SAAR 3896 segment 3, complete sequence', 'NC_038596': 'Lebombo virus isolate SAAR 3896 segment 4, complete sequence', 'NC_038597': 'Lebombo virus isolate SAAR 3896 segment 6, complete sequence', 'NC_038598': 'Lebombo virus isolate SAAR 3896 segment 8, complete sequence', 'NC_038599': 'Lebombo virus isolate SAAR 3896 segment 7, complete sequence', 'NC_038600': 'Lebombo virus isolate SAAR 3896 segment 2, complete sequence', 'NC_038601': 'Lebombo virus isolate SAAR 3896 segment 5, complete sequence', 'NC_038602': 'Lebombo virus isolate SAAR 3896 segment 10, complete sequence', 'NC_038603': 'Lebombo virus isolate SAAR 3896 segment 9, complete sequence', 'NC_038604': 'Orungo virus isolate UGMP 359 segment 1, complete sequence', 'NC_038605': 'Orungo virus isolate UGMP 359 segment 2, complete sequence', 'NC_038606': 'Orungo virus isolate UGMP 359 segment 4, complete sequence', 'NC_038607': 'Orungo virus isolate UGMP 359 segment 7, complete sequence', 'NC_038608': 'Orungo virus isolate UGMP 359 segment 8, complete sequence', 'NC_038609': 'Orungo virus isolate UGMP 359 segment 6, complete sequence', 'NC_038610': 'Orungo virus isolate UGMP 359 segment 3, complete sequence', 'NC_038611': 'Orungo virus isolate UGMP 359 segment 5, complete sequence', 'NC_038612': 'Orungo virus isolate UGMP 359 segment 9, complete sequence', 'NC_038613': 'Orungo virus isolate UGMP 359 segment 10, complete sequence', 'NC_038726': 'Catu virus strain BeH 151 nucleocapsid protein gene, complete cds', 'NC_038727': 'Catu virus strain BeH 151 polyprotein gene, complete cds', 'NC_038728': 'Catu virus strain BeH 151 RNA-dependent RNA polymerase gene, complete cds', 'NC_039024': 'Central cimpanzee simian foamy virus isolate BAD327, complete genome', 'NC_039025': 'Eastern chimpanzee simian foamy virus clone pHSRV13, complete genome', 'NC_039050': 'Cutavirus strain BR-337 NS1, putative VP1, hypothetical protein, VP2, and hypothetical protein genes, complete cds', 'NC_039191': 'Punta Toro virus strain Adames M protein gene, complete cds', 'NC_039192': 'Punta Toro virus strain Adames L protein gene, complete cds', 'NC_039193': 'Punta Toro virus strain Adames N protein and NS protein genes, complete cds', 'NC_039215': 'Cyclovirus VN isolate hcf2, complete genome', 'NC_036877': 'Cyclovirus PK5006, complete genome', 'NC_035889': 'Zika virus isolate ZIKV/H. sapiens/Brazil/Natal/2015, complete genome', 'NC_035758': 'Bastrovirus 7, complete genome', 'NC_035469': 'NY_014 poxvirus strain 2013, complete genome', 'NC_035475': 'Hudisavirus sp. isolate P22, complete genome', 'NC_034443': 'Le Dantec virus nucleoprotein, phosphoprotein, matrix, glycoprotein, hypothetical protein, and polymerase genes, complete cds', 'NC_034479': 'Bwamba virus strain M459 segment L polymerase gene, complete cds', 'NC_034480': 'Bwamba virus strain M459 segment S nucleoprotein and nonstructural protein NSs genes, complete cds', 'NC_034486': 'Guaroa virus strain BeH22063 segment S, complete sequence', 'NC_034487': 'Guaroa virus strain BeH22063 segment L, complete sequence', 'NC_034490': 'Bwamba virus strain M459 segment M glycoprotein precursor gene, complete cds', 'NC_034497': 'Madrid virus strain BT4075 RNA-dependent RNA polymerase gene, complete cds', 'NC_034498': 'Madrid virus strain BT4075 nucleocapsid protein and nonstructural protein genes, complete cds', 'NC_034505': 'Madrid virus strain BT4075 polyprotein gene, complete cds', 'NC_034506': 'Guaroa virus strain BeH22063 segment M, complete sequence', 'NC_034385': 'Human cosavirus F strain PK5006 polyprotein gene, complete cds', 'NC_034253': 'LI polyomavirus isolate LIPyV, complete genome', 'NC_032682': 'Indian encephalitis associated cyclovirus isolate IECSF08, complete genome', 'NC_032480': 'Husavirus sp. isolate 16370_59, complete genome', 'NC_030791': ' Hepatitis C virus genotype 7, complete genome', 'NC_030454': 'Enterovirus A114 strain V13-0285, partial genome', 'NC_030447': 'Gemycircularvirus HV-GcV1, complete genome', 'NC_030448': 'Gemycircularvirus HV-GcV2, complete genome', 'NC_030449': 'Unidentified circular ssDNA virus, complete genome', 'NC_030297': 'Torque teno mini virus 18 isolate 222, complete genome', 'NC_029647': 'Norovirus GIV, complete sequence', 'NC_026431': 'Influenza A virus (A/California/07/2009(H1N1)) segment 7 matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds', 'NC_026432': 'Influenza A virus (A/California/07/2009(H1N1)) segment 8 nuclear export protein (NEP) and nonstructural protein 1 (NS1) genes, complete cds', 'NC_026433': 'Influenza A virus (A/California/07/2009(H1N1)) segment 4 hemagglutinin (HA) gene, complete cds', 'NC_026434': 'Influenza A virus (A/California/07/2009(H1N1)) segment 6 neuraminidase (NA) gene, complete cds', 'NC_026435': 'Influenza A virus (A/California/07/2009(H1N1)) segment 2 polymerase PB1 (PB1) gene, complete cds; and nonfunctional PB1-F2 protein (PB1-F2) gene, complete sequence', 'NC_026436': 'Influenza A virus (A/California/07/2009(H1N1)) segment 5 nucleocapsid protein (NP) gene, complete cds', 'NC_026437': 'Influenza A virus (A/California/07/2009(H1N1)) segment 3 polymerase PA (PA) gene, complete cds', 'NC_026438': 'Influenza A virus (A/California/07/2009(H1N1)) segment 1 polymerase PB2 (PB2) gene, complete cds', 'NC_026422': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 1 polymerase PB2 (PB2) gene, complete cds', 'NC_026423': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 2 polymerase PB1 (PB1) and PB1-F2 protein (PB1-F2) genes, complete cds', 'NC_026424': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 3 polymerase PA (PA) and PA-X protein (PA-X) genes, complete cds', 'NC_026425': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 4 hemagglutinin (HA) gene, complete cds', 'NC_026426': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 5 nucleocapsid protein (NP) gene, complete cds', 'NC_026427': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 7 matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds', 'NC_026428': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 8 nuclear export protein (NEP) and nonstructural protein 1 (NS1) genes, complete cds', 'NC_026429': 'Influenza A virus (A/Shanghai/02/2013(H7N9)) segment 6 neuraminidase (NA) gene, complete cds', 'NC_025961': 'Cosavirus JMY-2014 isolate Cosa-CHN, complete genome', 'NC_025726': 'Torque teno mini virus ALA22, complete genome', 'NC_025727': 'Torque teno mini virus ALH8, complete genome', 'NC_025343': 'Sosuga virus isolate 2012, complete genome', 'NC_025114': 'Salivirus FHB, complete genome', 'NC_024888': 'Bufavirus-3 genes for NS1, putative VP1, hypothetical protein, VP2, complete cds, strain: BTN-63', 'NC_024781': 'Marburg marburgvirus isolate Ravn virus/H.sapiens-tc/KEN/1987/Kitum Cave-810040, complete genome', 'NC_024494': 'Heartland virus isolate Patient1 segment M, complete sequence', 'NC_024495': 'Heartland virus isolate Patient1 segment L, complete sequence', 'NC_024496': 'Heartland virus isolate Patient1 segment S, complete sequence', 'NC_024118': 'New Jersey polyomavirus-2013 isolate NJ-PyV-2013, complete genome', 'NC_024070': 'Rosavirus 2 strain GA7403, complete genome', 'NC_023888': 'Circo-like virus-Brazil hs1, complete genome', 'NC_022788': 'Gyrovirus Tu243, complete genome', 'NC_022789': 'Gyrovirus Tu789, complete genome', 'NC_022089': 'Parvovirus NIH-CQV putative 15-kDa protein, putative replication associated protein (rep), and putative capsid protein (cap) genes, complete cds', 'NC_020805': 'Chandipura virus isolate CIN 0451, complete genome', 'NC_020810': 'Duvenhage virus isolate 86132SA, complete genome', 'NC_020498': 'TTV-like mini virus isolate TTMV_LY1, complete genome', 'NC_020106': 'STL polyomavirus strain MA138, complete genome', 'NC_019843': 'Middle East respiratory syndrome-related coronavirus isolate HCoV-EMC/2012, complete genome', 'NC_019026': 'Astrovirus VA3 isolate VA3/human/Vellore/28054/2005, complete genome', 'NC_019027': 'Astrovirus VA4 isolate VA4/human/Nepal/s5363, complete genome', 'NC_019028': 'Astrovirus MLB3 isolate MLB3/human/Vellore/26564/2004, complete genome', 'NC_018401': 'Gyrovirus 4, complete genome', 'NC_018136': 'SFTS virus HB29 segment L, complete genome', 'NC_018137': 'SFTS virus HB29 segment S, complete genome', 'NC_018138': 'SFTS virus HB29 segment M, complete genome', 'NC_017091': 'Gyrovirus GyV3, complete genome', 'NC_018102': 'MW polyomavirus, complete genome', 'NC_016155': 'Astrovirus MLB2, complete genome', 'NC_015411': 'Sandfly Sicilian Turkey virus segment M, complete genome', 'NC_015412': 'Sandfly Sicilian Turkey virus segment L, complete genome', 'NC_015413': 'Sandfly Sicilian Turkey virus segment S, complete genome', 'NC_015373': 'Candiru virus segment M, complete genome', 'NC_015374': 'Candiru virus segment L, complete genome', 'NC_015375': 'Candiru virus segment S, complete genome', 'NC_014373': 'Bundibugyo ebolavirus, complete genome', 'NC_014361': 'Trichodysplasia spinulosa-associated polyomavirus, complete genome', 'NC_014372': 'Tai Forest ebolavirus isolate Tai Forest virus/H.sapiens-tc/CIV/1994/Pauleoula-CI, complete genome', 'NC_014089': 'Torque teno mini virus 5, complete genome', 'NC_014093': 'Torque teno midi virus 2, complete genome', 'NC_014095': 'Torque teno mini virus 6, complete genome', 'NC_014097': 'Torque teno mini virus 1, complete genome', 'NC_014068': 'Torque teno mini virus 8, complete genome', 'NC_013443': 'HMO Astrovirus A, complete genome', 'NC_013060': 'Astrovirus VA1, complete genome', 'NC_012986': 'Salivirus A isolate 02394-01, complete genome', 'NC_012957': 'Salivirus NG-J1, complete genome', 'NC_012798': 'Human cosavirus E1, complete genome', 'NC_012800': 'Cosavirus A strain HCoSV-A1 polyprotein gene, complete cds', 'NC_012802': 'Human cosavirus D1, complete genome', 'NC_012776': 'Lujo virus segment S, complete genome', 'NC_012777': 'Lujo virus segment L, complete genome', 'NC_009528': 'European bat lyssavirus 2 isolate RV1333, complete genome', 'NC_001798': 'Human alphaherpesvirus 2, complete genome', 'NC_027998': '67Human pegivirus 2, complete genome', 'NC_017993': 'Human papillomavirus 135, complete genome', 'NC_017994': 'Human papillomavirus 136, complete genome', 'NC_017995': 'Human papillomavirus type 137, complete genome', 'NC_017996': 'Human papillomavirus 140, complete genome', 'NC_017997': '71Human papillomavirus type 144, complete genome', 'NC_026946': 'Human papillomavirus KC5, complete genome', 'NC_021549': 'Human rotavirus B, complete genome', 'NC_027528': 'Human papillomavirus 201, complete genome', 'NC_021928': 'Human parainfluenza virus 4a, complete genome', 'NC_022518': '72Human endogenous retrovirus K113, complete genome', 'NC_007605': 'Human gammaherpesvirus 4, complete genome', 'NC_003461': 'Human respirovirus 1, complete genome', 'NC_009334': '764Human herpesvirus 4 type 2, complete genome', 'NC_006273': 'Human betaherpesvirus 5, complete genome', 'NC_027779': 'Human papillomavirus, complete genome', 'NC_001488': 'Human T-lymphotropic virus 2, complete genome', 'NC_028125': 'human papillomavirus 163, complete genome', 'NC_001405': 'Human mastadenovirus C, complete genome', 'NC_001722': '59Human immunodeficiency virus 2, complete genome', 'NC_002645': '7Human coronavirus 229E, complete genome', 'NC_005831': 'Human coronavirus NL63, complete genome', 'NC_000898': 'Human betaherpesvirus 6B, complete genome', 'NC_001716': 'Human betaherpesvirus 7, complete genome', 'NC_001796': 'Human respirovirus 3, complete genome', 'NC_001436': '7Human T-cell leukemia virus type I, complete genome', 'NC_011800': 'Human T-lymphotropic virus 4, complete genome', 'NC_001802': 'Human immunodeficiency virus 1, complete genome', 'NC_001806': 'Human alphaherpesvirus 1, complete genome', 'NC_003443': 'Human orthorubulavirus 2, complete genome', 'NC_001348': 'Human alphaherpesvirus 3, complete genome', 'NC_004104': 'Human papillomavirus type 90, complete genome', 'NC_038436': 'Human hepegivirus, complete genome', 'NC_039199': 'Human metapneumovirus, complete genome', 'NC_001664': '78Human betaherpesvirus 6A, complete genome', 'NC_009333': '69Human gammaherpesvirus 8, complete genome', 'NC_033781': 'Human papillomavirus type 156, complete genome', 'NC_006213': '741Human coronavirus OC43, complete genome', 'NC_006577': 'Human coronavirus HKU1, complete genome', 'NC_001781': 'Human orthopneumovirus, complete genome', 'NC_038235': 'Human orthopneumovirus, complete genome', 'NC_011202': 'Human mastadenovirus B, complete genome', 'NC_011203': 'Human mastadenovirus B, complete genome', 'NC_001460': 'Human mastadenovirus A, complete genome', 'NC_003266': 'Human mastadenovirus E, complete genome', 'NC_001454': 'Human mastadenovirus F, complete genome', 'NC_010956': 'Human mastadenovirus D, complete genome', 'NC_012959': 'Human adenovirus 54, complete genome', 'NC_035474': 'Human feces smacovirus 3, complete genome', 'NC_019023': 'human papillomavirus 166, complete genome', 'NC_023874': 'Human associated cyclovirus 10, complete genome', 'NC_023891': 'Human papillomavirus 178, complete genome', 'NC_023984': 'Human cosavirus, complete genome', 'NC_024472': 'Human astrovirus BF34, complete genome', 'NC_024694': 'Human circovirus VS6600022, complete genome', 'NC_020890': 'Human polyomavirus 12, complete genome', 'NC_021568': 'Human cyclovirus VS5700009, complete genome', 'NC_021542': 'Human rotavirus B, complete genome', 'NC_021483': '6Human papillomavirus 154, complete genome', 'NC_021541': 'Human rotavirus B, complete genome', 'NC_021543': '6Human rotavirus B, complete genome', 'NC_021544': 'Human rotavirus B, complete genome', 'NC_021545': '7Human rotavirus B, complete genome', 'NC_021546': 'Human rotavirus B, complete genome', 'NC_021547': 'Human rotavirus B, complete genome', 'NC_021550': 'Human rotavirus B, complete genome', 'NC_021551': 'Human rotavirus B, complete genome', 'NC_021548': '7Human rotavirus B, complete genome', 'NC_022095': 'Human papillomavirus 179, complete genome', 'NC_026252': 'Human smacovirus 1, complete genome', 'NC_026318': 'Human associated porprismacovirus 2, complete genome', 'NC_022892': 'Human papillomavirus 167, complete genome', 'NC_026817': 'Human genital-associated circular DNA virus-1, complete genome', 'NC_001355': 'Human papillomavirus type 6b, complete genome', 'NC_008188': 'human papillomavirus 103, complete genome', 'NC_004500': 'Human papillomavirus type 92, complete genome', 'NC_016157': '6Human papillomavirus 126, complete genome', 'NC_010329': '6Human papillomavirus type 88, complete genome', 'NC_007018': 'Human parvovirus 4 G1, complete genome', 'NC_001593': '56Human papillomavirus type 53, complete genome', 'NC_012485': 'Human papillomavirus 109, complete genome', 'NC_012486': '7Human papillomavirus type 112, complete genome', 'NC_007026': 'Human picobirnavirus, segment 1, complete sequence', 'NC_007027': 'Human picobirnavirus, segment 2, complete sequence', 'NC_000883': 'Human parvovirus B19, complete genome', 'NC_001591': 'Human papillomavirus type 49, complete genome', 'NC_001458': 'Human papillomavirus type 63, complete genome', 'NC_001595': '7Human papillomavirus type 7, complete genome', 'NC_001531': 'Human papillomavirus 5, complete genome', 'NC_013035': 'Human papillomavirus 116, complete genome', 'NC_001526': '6Human papillomavirus type 16, complete genome', 'NC_001538': 'Human polyomavirus 1, complete genome', 'NC_014955': 'Human papillomavirus 132, complete genome', 'NC_001596': 'Human papillomavirus 9, complete genome', 'NC_001576': '9Human papillomavirus type 10, complete genome', 'NC_001583': '55Human papillomavirus type 26, complete genome', 'NC_001586': '61Human papillomavirus type 32, complete genome', 'NC_001587': 'Human papillomavirus type 34, complete genome', 'NC_014185': 'Human papillomavirus 121, complete genome', 'NC_010810': '61Human TMEV-like cardiovirus, complete genome', 'NC_014406': 'Human polyomavirus 6, complete genome', 'NC_014407': 'Human polyomavirus 7, complete genome', 'NC_012729': 'Human bocavirus 4 NI, complete genome', 'NC_008189': 'Human papillomavirus type 101, complete genome', 'NC_009238': '0Human polyomavirus 3, complete genome', 'NC_012564': 'Human bocavirus 3, complete genome', 'NC_014469': 'Human papillomavirus 127, complete genome', 'NC_012042': 'Human bocavirus 2c PK, complete genome', 'NC_012801': 'Human cosavirus B, complete genome', 'NC_001690': 'Human papillomavirus type 48, complete genome', 'NC_001693': 'Human papillomavirus type 60, complete genome', 'NC_001691': 'Human papillomavirus type 50, complete genome', 'NC_001676': 'Human papillomavirus type 54, complete genome', 'NC_005134': 'Human papillomavirus type 96, complete genome', 'NC_004295': 'Human erythrovirus V9, complete genome', 'NC_001457': 'Human papillomavirus 4, complete genome', 'NC_015150': '6Human polyomavirus 9, complete genome', 'NC_001354': 'Human papillomavirus type 41, complete genome', 'NC_001943': 'Human astrovirus, complete genome', 'NC_038392': 'Human stool-associated circular virus NG13, complete genome', 'NC_038497': 'Human gemycircularvirus GeTz1, complete genome', 'NC_038522': 'Human papillomavirus type 161, complete genome', 'NC_038523': 'Human papillomavirus 172, complete genome', 'NC_038524': 'Human papillomavirus 175, complete genome', 'NC_038525': '7Human papillomavirus 204, complete genome', 'NC_038878': 'Human rhinovirus NAT001, complete genome', 'NC_038889': '52human papillomavirus 30, complete genome', 'NC_038914': 'Human papillomavirus 184, complete genome', 'NC_039061': 'Human associated huchismacovirus 1, complete genome', 'NC_039062': 'Human associated huchismacovirus 2, complete genome', 'NC_039063': 'Human associated huchismacovirus 3, complete genome', 'NC_039070': 'Human feces smacovirus 2, complete genome', 'NC_039086': 'Human papillomavirus 187, complete genome', 'NC_039089': 'human papillomavirus 71, complete genome', 'NC_038319': 'Human parechovirus 1, complete genome', 'NC_034616': 'Human papillomavirus type 85, complete genome', 'NC_028459': 'Human associated gemyvongvirus 1, complete genome', 'NC_012213': 'human papillomavirus 108, complete genome', 'NC_014952': 'Human papillomavirus type 128, complete genome', 'NC_014953': 'Human papillomavirus type 129, complete genome', 'NC_014954': 'Human papillomavirus type 131, complete genome', 'NC_014956': 'Human papillomavirus 134, complete genome', 'NC_035211': 'Human fecal virus Jorvi4, complete genome', 'NC_035213': 'Human fecal virus Jorvi3, complete genome', 'NC_035212': 'Human fecal virus Jorvi2, complete genome', 'NC_040640': '72Human papillomavirus type 203, complete genome', 'NC_040306': 'Human feces pecovirus, complete genome', 'NC_040309': '5Human DNA virus, complete genome', 'NC_055523': '6Human lung-associated vientovirus FB, complete genome', 'NC_038882': ' Hepatitis C virus (isolate H77) genotype 1, complete cds', 'NC_004102': ' Hepatitis C virus genotype 1, complete genome', 'NC_009824': ' Hepatitis C virus genotype 3, genome', 'NC_009827': ' Hepatitis C virus genotype 6, complete genome', 'NC_009823': ' Hepatitis C virus genotype 2, complete genome', 'NC_009826': ' Hepatitis C virus genotype 5, genome', 'NC_009825': ' Hepatitis C virus genotype 4, genome', 'NC_003977': ' Hepatitis B virus (strain ayw) genome', 'AC_000007': ' Human adenovirus 2, complete genome', 'AC_000008': ' Human adenovirus 5, complete genome', 'AC_000019': ' Human adenovirus type 35, complete genome', 'AC_000017': ' Human adenovirus type 1, complete genome', 'AC_000018': ' Human adenovirus type 7, complete genome', 'AC_000006': ' Human adenovirus D, complete genome', 'AC_000005': ' Human mastadenovirus A, complete genome', 'NC_001401': 'Adeno-associated virus - 2, complete genome', 'J04353': ' Human papillomavirus type 31 (HPV-31) complete genome', 'X74481': ' Human papillomavirus type 52 genomic DNA', 'X74479': ' Human papillomavirus type 45 genomic DNA', 'DQ080079 ': ' Human papillomavirus type 68a, complete genome', 'X05015': ' Human papillomavirus type 18 E6, E7, E1, E2, E4, E5, L1 & L2 genes', 'M12732': ' Human papillomavirus type 33, complete genome', 'X94165': 'Human papillomavirus type 73 E6, E7, E1, E2, E4, L2, and L1 genes', 'X74477': ' Human papillomavirus type 35H genomic DNA', 'V01116': 'Human papilloma virus type 1a E6, E8, E7, L2 & L2 genes', 'U21941': 'Human papillomavirus type 70, complete genome', 'X74472': 'Human papillomavirus type 26 genomic DNA', 'X77858': 'Human papilloma virus type 59, complete viral genome', 'M62849': 'Human papillomavirus ORFs', 'M62877': 'Human papilloma virus type 51 genomic DNA, partial sequence', 'X74474': ' Human papillomavirus type 30 genomic DNA', 'D90400': 'Human papillomavirus type 58 complete genome'}


    f2 = open(infile, 'r')
    l = f2.readlines()
    names = [l[i].split('\n')[0].split(':')[0].replace('>', '') for i in range(1, len(l) - 2, 2)]
    seqs = [handle_non_ATGC(l[i].split('\n')[0]) for i in range(2, len(l) - 1, 2)]

    mapper = []
    for i in range(len(fasta_id)):
        id = findi(names, fasta_id[i])
        mapper.append(id[0])

    all_len = [len(seqs[index]) for index in mapper]
    all_con = [names[index] for index in mapper]

    dic_con= dict(zip(all_con, all_len))


    contig = []
    seq_len = []
    average = []
    cover_av = []
    cover_max = []
    full_vir2 = []
    for i in range(len(uni_vir)):
        indexes = findi(vir, uni_vir[i])

        con = list(set([all_con[k] for k in indexes]))
        avg = np.mean([(score[k]) for k in indexes])
        average.append(avg)
        cover_av.append(np.mean([(coverage[k]) for k in indexes]))
        cover_max.append(max([(coverage[k]) for k in indexes]))
        contig.append(con)
        seq_len.append(get_value(dic_con, con))
        if uni_vir[i] in list(ref_dict.keys()):
            full_vir2.append(ref_dict[uni_vir[i]])
        else:
            full_vir2.append('other')

    smp_name = [infile.split('/')[-1].split('.')[0].replace('_gdc_realn_rehead_contigs', '') for i in uni_vir]

    case_id, sample_tp, project_id = get_project_id(smp_name)

    df = pd.DataFrame({'case_id': case_id, 'sample': smp_name, 'virus': uni_vir,'full_virus2': full_vir2,
                       'avg_similarity': average, 'avg_coverage': cover_av, 'max_coverage': cover_max, 'lengths': seq_len, 'contig': contig,
                      'project_id': project_id, 'sample_type': sample_tp})

    df = df.set_index('case_id')

    idx = np.unique([i for i in range(len(df.index)) if (list(df['virus'])[i] not in contaminators) and (list(df['max_coverage'])[i] >= 90) and (list(df['avg_coverage'])[i] > 50) and (list(df['avg_similarity'])[i] > 90)])

    df = df.iloc[idx]

    df.to_csv(filen + '.csv')

    return df

def get_table_novel(infile,outfie,filen,nuc_len_translate = 100):
    '''
    :param infile: input file (contigs-fasta), at dbgap-mutporg/
    :param outfie: An output file from input file
    :param filen: generated files name (dirname_ref+'/'+file.split('/')[-1].split('.')[0])
    :return: csv data with refseq families
   '''
    findi = lambda x, v: [i for i, j in enumerate(x) if j == v]
    flatten = lambda l: [item for sublist in l for item in sublist]
    get_value = lambda x, y: [x.get(i) for i in y]

    col_names = ["id", "virus", "score", "coverage","b","c","d","e","f","g","h","i"]
    out11=pd.read_csv(outfie, sep='\t', header=None, names=col_names)

    fasta_ids = list(out11['id'])
    fasta_id = [i.split('_bact')[0] for i in fasta_ids]

    virs = list(out11['virus'])
    vir=[i.split('|')[2].split('.')[0] for i in virs]

    uni_vir = sorted(list(set(vir)))
    score = list(out11['score'])
    coverage = list(out11['coverage'])

    with open('virus_names.json') as fp:
        virus_dict = json.load(fp)

    f2 = open(infile,'r')
    l = f2.readlines()
    names=[l[i].split('\n')[0].split('_bact')[0].replace('>','') for i in range(1,len(l)-2,2)]
    seqs = [handle_non_ATGC(l[i].split('\n')[0]) for i in range(2, len(l)-1, 2)]

    mapper=[]
    for i in range(len(fasta_id)):
        id = findi(names, fasta_id[i])

        mapper.append(id[0])

    all_len = [len(seqs[index]) for index in mapper]
    all_con = [names[index] for index in mapper]
    all_seq = [seqs[index] for index in mapper]

    dic_con = dict(zip(all_con, all_len))

    max_prot = []
    for n in range(len(all_seq)):
        beans = all_seq[n]
        if len(all_seq[n]) > nuc_len_translate:
            x1 = translate_frameshifted(beans[0:])
            x2 = translate_frameshifted(beans[1:])
            x3 = translate_frameshifted(beans[2:])
            x4 = translate_frameshifted(reverse_complement(beans))
            x5 = translate_frameshifted(reverse_complement(beans[:len(beans) - 1]))
            x6 = translate_frameshifted(reverse_complement(beans[:len(beans) - 2]))

            x = [x1, x2, x3, x4, x5, x6]
            l = flatten([i.split('_') for i in x])
            max_prot.append(max(len(n) for n in l))

        else:
            max_prot.append(0)

    average=[]
    cover_av = []
    full_vir=[]
    contig = []
    max_protein = []
    seq_len = []
    for i in range(len(uni_vir)):
        indexes = findi(vir, uni_vir[i])
        avg=np.mean([(score[k]) for k in indexes])
        con = list(set([all_con[k] for k in indexes]))
        m_p = [max_prot[k] for k in indexes]
        average.append(avg)
        cover_av.append(np.mean([(coverage[k]) for k in indexes]))
        full_vir.append(virus_dict[uni_vir[i]])
        contig.append(con)
        seq_len.append(get_value(dic_con, con))
        max_protein.append(m_p)


    smp_name = [infile.split('/')[-1].split('.')[0].replace('_gdc_realn_rehead_contigs','') for i in uni_vir]

    case_id, sample_tp, project_id = get_project_id(smp_name)

    df=pd.DataFrame({'case_id':case_id, 'sample':smp_name,'virus':uni_vir,'full_virus':full_vir,
                     'avg_similarity':average, 'avg_coverage': cover_av,'lengths': seq_len, 'contig': contig,
                     'max_protein': max_protein,'project_id':project_id,'sample_type':sample_tp})


    df = df.set_index('case_id')
    idx = np.unique([i for i in range(len(df.index)) if
                     (((np.max(df['lengths'][i])>300 and list(df['avg_coverage'])[i] > 40)  or (list(df['avg_coverage'])[i] > 75) and
                                                           list(df['avg_similarity'])[i] > 75)) and (list(df['virus'])[i] not in contaminators) and np.max(df['lengths'][i])>100])

    df = df.iloc[idx]

    df.to_csv(filen+'.csv')

    return df


def combine_csv(indir,outname,extension = 'csv'):
    '''
    combine all csv files to one file and remove them
    :param indir:directory (dirname_ref or dirname_herv)
    :param outname: name of the combined csv file
    :return: combined csv file
    '''
    all_filenames = [i for i in glob.glob(indir + '/*.{}'.format(extension))]
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    frm = glob.glob(indir +'/' + '*.csv')
    [os.system('rm -rf ' + frm[i]) for i in range(len(frm))]

    combined_csv.to_csv(outname, index=False, encoding='utf-8-sig')


def get_blast_results(params):
    '''

    :param output_dir:
    :return:
    '''
    output_dir = params[0]
    file = params[1]

    fn = dirname_novel+'/'+file.split('/')[-1].split('.')[0] +'.csv'
    if not exists(fn):
        subprocess.run('blastn -db $BLASTDB/C-RVDBv16.0.txt -query '+file+'  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue" -max_target_seqs 1 -evalue 0.01 | sort -u > '+output_dir+'tmpout_novel'+file.split('/')[-1].split('.')[0]+'.txt',shell=True)
        get_table_novel(file,output_dir+'tmpout_novel'+file.split('/')[-1].split('.')[0]+'.txt',dirname_novel+'/'+file.split('/')[-1].split('.')[0])
        subprocess.run('rm -rf '+output_dir+'tmpout_novel'+file.split('/')[-1].split('.')[0]+'.txt',shell=True)

    fn = dirname_ref + '/' + file.split('/')[-1].split('.')[0]
    if not exists(fn):
        subprocess.run('blastn -db $BLASTDB/human_virus_db.fa -query ' + file + '  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue"  -evalue 0.01  -word_size 15 | sort -u >  '+output_dir+'tmpout_refseq'+file.split('/')[-1].split('.')[0]+'.txt',shell=True)
        get_table_reference(file, output_dir+'tmpout_refseq'+file.split('/')[-1].split('.')[0]+'.txt', dirname_ref + '/' + file.split('/')[-1].split('.')[0])
        subprocess.run('rm -rf '+output_dir+'tmpout_refseq'+file.split('/')[-1].split('.')[0]+'.txt',shell=True)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--input", type=str, help="input directory", required=True
    )
    parser.add_argument(
        "--output", type=str, help="output directory", required=False
    )
    args = parser.parse_args()

    if args.output is None:
        output_dir = args.input
        output_dir = (os.getcwd()+'/'+'outputs'+'/')

    else:
        output_dir = args.output+'outputs'+'/'

    os.system('mkdir ' + output_dir[:-1])

    os.system('mkdir ' + output_dir + 'out_novel')
    os.system('mkdir ' + output_dir + 'out_ref')

    dirname_ref = output_dir + 'out_ref'
    dirname_novel = output_dir + 'out_novel'

    files = glob.glob(args.input + '*.txt')

    freeze_support()
    pool = mp.Pool(processes=20)
    pool.map(get_blast_results, [[output_dir, f] for f in files if f is not None])

    combine_csv(dirname_ref, dirname_ref + "/combined_csv_ref.csv")
    combine_csv(dirname_novel, dirname_novel + "/combined_csv_novel.csv")