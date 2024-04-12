<!DOCTYPE html>
<html>

<head>
	<style>
	table, th, td {
    	border: 1px solid black;
    	border-collapse: collapse;
	}
	th, td {
    	padding: 5px;
	}
	</style>
</head>

<body>


<?php

$w = stream_get_wrappers();
echo 'openssl: ',  extension_loaded  ('openssl') ? 'yes':'no', "\n";
echo 'http wrapper: ', in_array('http', $w) ? 'yes':'no', "\n";
echo 'https wrapper: ', in_array('https', $w) ? 'yes':'no', "\n";
echo 'wrappers: ', var_export($w);

//Define global variables:
$singleAA = array('G','P','A','V','L','I','M','F','Y','W','S','T','C','N','Q','D','E','K','R','H');
$molwAA = array(
	'G'=>75.0669,
	'P'=>115.1310,
	'A'=>89.0935,
	'V'=>117.1469,
	'L'=>131.1736,
	'I'=>131.1736,
	'M'=>149.2124,
	'F'=>165.1900,
	'Y'=>181.1894,
	'W'=>204.2262,
	'S'=>105.0930,
	'T'=>119.1197,
	'C'=>121.1590,
	'N'=>132.1184,
	'Q'=>146.1451,
	'D'=>133.1032,
	'E'=>147.1299,
	'K'=>146.1882,
	'R'=>174.2017,
	'H'=>155.1552
	);
	//http://www.webqc.org/aminoacids.php

$proase = $_GET["proase"];
$url = $_GET["url"];
$subpro = $_GET["subpro"];
$sequence = str_split(strtoupper($_GET["sequence"]));
// $startp = $_GET["startp"];
// $endp = $_GET["endp"];

// $proase = "catS";
// $url = "C:\Users\davisjad\OneDrive - National Institutes of Health\merops\merops.html";
// $subpro = "catK";
// $sequence = "AAAAAAAA";
$startp = 1;
$endp = 8;

$title = $proase." on ".$subpro;

$scrape = file_get_contents($url);
$tablestr = stristr($scrape,'summary="Specificity matrix"');
$tablestr = strstr($tablestr,'</table>',TRUE);

$isAll = 'yay';
foreach ($sequence as $possAA) {
	if (count(array_search($possAA,$singleAA))==0) {
		$isAll = 'nay';
	}
}

if (empty($tablestr)) {
	print("Protease-ase link does not contain a specificity matrix.");
} elseif ($isAll=='nay') {
	print("Substrate sequence is not valid. Only use single-letter amino acid abbreviations.");
} elseif ($endp-$startp<=0) {
	print("Upper bound must be larger than lower bound.");
} else {
	$scrapeAA = array('Gly','Pro','Ala','Val','Leu','Ile','Met','Phe','Tyr','Trp','Ser','Thr','Cys','Asn','Gln','Asp','Glu','Lys','Arg','His');
	$numAA = count($scrapeAA);

	$arr = array_fill(0,$numAA,0);

	for ($n=0;$n<$numAA;$n++) {
		$lookAA = $scrapeAA[$n];
		//Isolate the row for the current amino acid:
		$currstr = strstr($tablestr,$lookAA);
		$currstr = strstr($currstr,"<td");
		$currstr = strstr($currstr,"</tr>",TRUE);

		$currstr = str_replace("</td>",";",$currstr);
		$currstr = strip_tags($currstr);
		$currstr = str_replace("\n","",$currstr);
		$currstr = str_replace(" ","",$currstr);

		$currarr = [];
		$token = strtok($currstr,";");

		while ($token !== false) {
			$currarr = array_merge($currarr,[(float)$token]);
			$token = strtok(";");
		}

		$arr[$n] = $currarr;
	}

	// /*
	// Dummy variables:
	$arr = array(
		array(3,3,0,0,0,0,0,4),
		array(1,2,0,0,0,0,1,0),
		array(0,0,0,2,0,0,1,1),
		array(1,1,2,1,0,0,0,1),
		array(0,1,6,1,1,1,0,1),
		array(0,1,0,0,0,0,0,0),
		array(0,0,0,1,0,0,0,0),
		array(0,1,4,1,0,1,1,0),
		array(0,0,0,0,0,0,4,0),
		array(0,0,0,0,0,0,0,0),
		array(1,0,0,1,0,0,0,0),
		array(0,0,0,0,1,0,0,0),
		array(0,0,0,0,0,0,0,0),
		array(0,0,0,0,0,0,0,0),
		array(0,0,0,1,0,0,1,0),
		array(0,0,0,0,0,0,0,1),
		array(1,1,0,0,1,1,0,1),
		array(1,0,1,1,5,1,1,0),
		array(0,0,1,5,1,5,0,0),
		array(1,0,0,0,0,0,0,0)
		);
	$sequence = array('V','P','M','S','M','R','G','G');
	$startp = 3;
	$endp = 6;
	// */

	analyze($arr, $sequence, $startp, $endp);
}


function AA2molw($AA) {
	return $GLOBALS['molwAA'][$AA];
}

function analyze($arr, $sequence, $startp, $endp) {
	$seqlen = count($sequence);
	$testlen = $endp-$startp+1;
	$numseq = $seqlen-$testlen+1;
	$listAA = $GLOBALS['singleAA'];

	//Convert sequence of AAs to an array of molecular weights:
	$seqmolw = array_map('AA2molw',$sequence);

	$fircol = array_column($arr,0);
	$fircolsum = array_sum($fircol);
	$seccol = array_column($arr,1);
	$seccolsum = array_sum($seccol);
	$totNumSubst = array($fircolsum,$seccolsum);
	
	for ($i=2;$i<count($arr[0]);$i++) {
		$col = array_column($arr,$i);
		$colsum = array_sum($col);
		$totNumSubst = array_merge($totNumSubst,array($colsum));
	}

	$scores = array_fill(0,$numseq,0);
	$seqArr = array_fill(0,$numseq,0);
	$scoArr = array_fill(0,$numseq,0);
	$sNoArr = array_fill(0,$numseq,0);
	$sPoArr = array_fill(0,$numseq,0);
	$ePoArr = array_fill(0,$numseq,0);
	$eScArr = array_fill(0,$numseq,0);
	$rScArr = array_fill(0,$numseq,0);
	$mScArr = array_fill(0,$numseq,0);
	$lCLArr = array_fill(0,$numseq,0);
	$rCLArr = array_fill(0,$numseq,0);
	$lCWArr = array_fill(0,$numseq,0);
	$rCWArr = array_fill(0,$numseq,0);


	for ($i=0;$i<$numseq;$i++) {
		$curseq = array_slice($sequence,$i,$testlen);
		$curscore = 0;
		$curseqscore = array_fill(0,$testlen,0);
		$curseqPscore = array_fill(0,$testlen,0);

		for ($j=0;$j<$testlen;$j++) {
			$whichAA = array_search($curseq[$j],$listAA);
			$correctedJ = $startp+$j-1;
			$curscore = $curscore+$arr[$whichAA][$correctedJ];
			$curseqscore[$j] = $arr[$whichAA][$correctedJ];
			$curseqPscore[$j] = $arr[$whichAA][$correctedJ]/$totNumSubst[$correctedJ];
		}

		$sortedScore = $curseqscore;
		sort($sortedScore);
		$minScore = $sortedScore[0];
		$maxScore = $sortedScore[$testlen-1];
		$range = [$minScore,$maxScore];
		if (($testlen/2-floor($testlen/2))!=0) {
			$medScore = $sortedScore[floor($testlen/2)];
		} else {
			$medScore = ($sortedScore[floor($testlen/2)-1]+$sortedScore[floor($testlen/2)])/2;
		}
		$Snormi = array_sum($curseqPscore);

		$lclen = $i+($testlen/2);
		$rclen = $seqlen-($i+($testlen/2));
		$lcweight = array_slice($seqmolw,0,$lclen);
		$rcweight = array_slice($seqmolw,$lclen,$rclen);

		//Make basically a structure array:
		$scores[$i] = array(
			'sequence'=>$curseq,
			'score'=>$curscore,
			'scoreNorm'=>$Snormi,
			'startpos'=>$i+1,
			'endpos'=>$i+$testlen,
			'eachScore'=>$curseqscore,
			'rangeScore'=>$range,
			'medianScore'=>$medScore,
			'leftcleavagelen'=>$lclen,
			'rightcleavagelen'=>$rclen,
			'leftcleavageweight'=>(array_sum($lcweight)-($lclen-1)*18.015)/1000,
			'rightcleavageweight'=>(array_sum($rcweight)-($rclen-1)*18.015)/1000
			);

		$seqArr[$i] = $curseq;
		$scoArr[$i] = $curscore;
		$sNoArr[$i] = $Snormi;
		$sPoArr[$i] = $i+1;
		$ePoArr[$i] = $i+$testlen;
		$eScArr[$i] = $curseqscore;
		$rScArr[$i] = $range;
		$mScArr[$i] = $medScore;
		$lCLArr[$i] = $lclen;
		$rCLArr[$i] = $rclen;
		$lCWArr[$i] = (array_sum($lcweight)-($lclen-1)*18.015)/1000;
		$rCWArr[$i] = (array_sum($rcweight)-($rclen-1)*18.015)/1000;
	}

	array_multisort($scoArr,SORT_DESC,$seqArr,$sNoArr,$sPoArr,$ePoArr,$eScArr,$rScArr,$mScArr,$lCLArr,$rCLArr,$lCWArr,$rCWArr);

	//Make HTML code:
	$HTML =
	"<table style=\"width:100%\">\n
	<tr>\n
	<th colspan=\"13\">".$GLOBALS['title']."</th>\n
	</tr>\n
	<tr>\n
	<th>Rank</th>\n
	<th>Cumulative Score</th>\n
	<th>Normal Score</th>\n
	<th>Start</th>\n
	<th>Sequence</th>\n
	<th>End</th>\n
	<th>Scores</th>\n
	<th>Range</th>\n
	<th>Mean</th>\n
	<th>Left Length (A.A.)</th>\n
	<th>Right Length (A.A.)</th>\n
	<th>Left Weight (kDa) </th>\n
	<th>Right Weight (kDa) </th>\n
	</tr>\n\n";
	/*
	foreach ($scores as $currow) {
		$currHTML =
		"<tr>\n
		<td>".(string)$currow['score']."</td>\n
		<td>".(string)round($currow['scoreNorm'],3)."</td>\n
		<td>".(string)$currow['startpos']."</td>\n
		<td>".implode("",$currow['sequence'])."</td>\n
		<td>".(string)$currow['endpos']."</td>\n
		<td>"."[ ".implode(" ",$currow['eachScore'])." ]"."</td>\n
		<td>"."[ ".implode(" ",$currow['rangeScore'])." ]"."</td>\n
		<td>".(string)$currow['medianScore']."</td>\n
		<td>".(string)$currow['leftcleavagelen']."</td>\n
		<td>".(string)$currow['rightcleavagelen']."</td>\n
		<td>".(string)round($currow['leftcleavageweight'],3)."</td>\n
		<td>".(string)round($currow['rightcleavageweight'],3)."</td>\n";
		$HTML = $HTML.$currHTML."</tr>\n\n";
	}
	*/
	for ($k=0;$k<$numseq;$k++) {
		$currHTML =
		"<tr>\n
		<td>".(string)($k+1)."</td>\n
		<td>".(string)$scoArr[$k]."</td>\n
		<td>".(string)round($sNoArr[$k],3)."</td>\n
		<td>".(string)$sPoArr[$k]."</td>\n
		<td>".implode("",$seqArr[$k])."</td>\n
		<td>".(string)$ePoArr[$k]."</td>\n
		<td>"."[ ".implode(" ",$eScArr[$k])." ]"."</td>\n
		<td>"."[ ".implode(" ",$rScArr[$k])." ]"."</td>\n
		<td>".(string)$mScArr[$k]."</td>\n
		<td>".(string)$lCLArr[$k]."</td>\n
		<td>".(string)$rCLArr[$k]."</td>\n
		<td>".(string)round($lCWArr[$k],3)."</td>\n
		<td>".(string)round($rCWArr[$k],3)."</td>\n";
		$HTML = $HTML.$currHTML."</tr>\n\n";
	}
	$HTML = $HTML."</table>\n";

	print($HTML);
}

?>


</body>

</html>
