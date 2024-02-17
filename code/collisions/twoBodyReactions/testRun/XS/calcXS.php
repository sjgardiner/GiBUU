<?php

//    echo '<pre>';
//    print_r($_POST);
//    echo '</pre>';

//Array
//(
//    [ID1] => 1
//    [anti1] => false
//    [ID2] => 1
//    [anti2] => false
//    [Charge1] => 0
//    [Charge2] => 0
//)

$path = getcwd();
chdir("run");

$output = shell_exec('/usr/bin/rm -f *dat');

//===== generate the jobcard =====

$file = fopen("jobCard","w");
fwrite($file,"\$Plotter".PHP_EOL);
fwrite($file,"id1=".$_POST['ID1'].PHP_EOL);
fwrite($file,"id2=".$_POST['ID2'].PHP_EOL);
fwrite($file,"q1=".$_POST['Charge1'].PHP_EOL);
fwrite($file,"q2=".$_POST['Charge2'].PHP_EOL);
fwrite($file,"anti1=".$_POST['anti1'].PHP_EOL);
fwrite($file,"anti2=".$_POST['anti2'].PHP_EOL);

$file_data = file_get_contents('../jobCard.tail');
fwrite($file,$file_data.PHP_EOL);

fclose($file);

//===== run GiBUU =====


//echo "running GiBUU...";
$output = shell_exec('../../CrossSectionPlotter.x < jobCard > run.log 2>&1');
//echo "<pre>$output</pre>";
//echo "finished?";

//$output = shell_exec('ls -lart');
//echo "<pre>$output</pre>";

$file_data = file_get_contents('XS.html');
echo $file_data;

//===== run gnuplot =====

$output = shell_exec('/usr/bin/rm -f plot*.svg');

if (file_exists('XS.dat'))
{
    $output = shell_exec('gnuplot ../plot.gp > gnuplot.log 2>&1');
    //echo "<pre>$output</pre>";
    //echo "finished?".PHP_EOL;
}

//===== embed plots =====

if (file_exists('plot.svg'))
{
    echo file_get_contents("plot.svg")."<br>".PHP_EOL;
    echo file_get_contents("plot_plab.svg")."<br>".PHP_EOL;

    //    echo "<img src=\"plot.svg\" type=\"image/svg+xml\" alt=\"Plot\">".PHP_EOL;
    //    echo "<embed src=\"file://$path/run/plot.svg\" type=\"image/svg+xml\" alt=\"Plot\"></embed>".PHP_EOL;
    //    echo "<embed src=\"file://$path/run/plot_plab.svg\" type=\"image/svg+xml\" alt=\"Plot\"></embed>".PHP_EOL;
}
else
{
    echo "<h2>Sorry, wrong input ...</h2>".PHP_EOL;
}

$file_data = file_get_contents('../foot.html');
echo $file_data;

//===== finish =====

chdir($path);



?>
