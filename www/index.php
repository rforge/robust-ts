
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<hr style="width: 100%; height: 4px;">
<h2>Guidelines for package <span
style="font-family: Courier New,Courier,monospace;">robust-ts</span>:&nbsp;</h2>
<h3>naming convention&nbsp;</h3>
<ul>
<li><span style="font-weight: normal;">yet to be announced (to be
decided by plenary vote in r-sig-robust)<br>
</span></li>
<li><span style="font-weight: normal;">open issues</span></li>
<ul>
<li><span style="font-weight: normal;">package &amp; function
names: suffix vs. prefix ``rob'' / </span><span
style="font-weight: normal;">``robust''</span></li>
<li><span style="font-weight: normal;">captialization, separators<br>
</span></li>
</ul>
<li><span style="font-weight: normal;">part <span
style="font-family: Courier New,Courier,monospace;">ts</span> is
agreed upon<br>
</span></li>
</ul>
<h3>programming in layers/stages<br>
</h3>
<ul>
<li>
<h4>first layer/stage (implementation of algorithms)<br>
</h4>
</li>
<ul>
<li>implementation of the algorithm in R</li>
<li>if necessary: delegation of critical code to C / C++ or FORTRAN</li>
<li>basic interfacing functions (with large number of arguments)<br>
</li>
<li>no particular data types (apart from vectors and&nbsp; matrices)</li>
<li>common "style" not necessary<br>
</li>
<li>no encapsulation</li>
<li>gives "alpha verion" --- extensive testing at this stage
possible and necessary<br>
</li>
</ul>
<li>
<h4>second layer/stage (user interface)</h4>
</li>
<ul>
<li>for consistency within this package: common "style" necessary </li>
<li>use of particular classes for input and output<br>
</li>
<li>low number of structured arguments in user-interface
functions/methods<br>
</li>
<li>mimic existing user interfaces</li>
<ul>
<li>use of data_frames in regression type setup<br>
</li>
</ul>
<ul>
<li>use of formula type model definitions</li>
</ul>
<li>accessor and replacment functions for the slots of class<br>
</li>
<li>particular methods <br>
</li>
<ul>
<li>print <br>
</li>
<li>show</li>
<li>plot</li>
<li>summary<br>
</li>
</ul>
</ul>
</ul>
<h3>input structure </h3>
<ul>
<li> possible input data: <br>
</li>
<ul>
<li>(raw) <span style="font-family: Courier New,Courier,monospace;">vector</span>s,
<br>
</li>
<li><span style="font-family: Courier New,Courier,monospace;">ts</span>
objects</li>
<li><span style="font-family: Courier New,Courier,monospace;">zoo</span>
objects (if method allows for irregularly spaced time stamps)<br>
</li>
</ul>
<li> function argument list / signature: <br>
as far as possible same arg's as classical routine </li>
<li> additional arguments split into control/``design
type''/method type arguments:
<ul>
<li> ``design'':&nbsp; tuning parameters for the robust procedure
which are <br>
</li>
<ul>
<li>subject to frequent change by the standard
user; <br>
</li>
<li>explicit and transparent to the user <br>
</li>
<li>e.g. window width in some filters<br>
</li>
</ul>
<li> <span style="font-family: Courier New,Courier,monospace;">control</span>:
tuning parameters for the robust procedure&nbsp;&nbsp;</li>
<ul>
<li>for which the standard user will use default values<br>
</li>
<li>these are packed into control object</li>
<li>may be generated by a function <span
style="font-family: Courier New,Courier,monospace;">control</span>
--see <span style="font-family: Courier New,Courier,monospace;">robustbase</span>
</li>
<li>the author will provide control-functions for each of
his provided methods (borrow from <span
style="font-family: Courier New,Courier,monospace;">robustbase</span>)</li>
</ul>
<li> <span style="font-family: Courier New,Courier,monospace;">method</span>
argument: used in a ``super''-function to determine the actual
robustification / procedure<br>
two different types of values should be allowed simultaneously
<ul>
<li> (a vector of) <span
style="font-family: Courier New,Courier,monospace;">character</span>s:
<br>
</li>
<ul>
<li>for: default/standard
user</li>
<li>names of the robustification(s)/ procedure(s) </li>
<li>vector-valued: several methods may be calculated in
parallel / in one
loop
over the time stamps </li>
</ul>
<li>an object of a S4-method-class providing the functions to
be used for this/these actual procedure(s)<br>
</li>
<ul>
<li>for: ``power
user''<br>
</li>
<li>goal:
a ``power user'' should be able to define his own method-class<br>
to make the routine use the user's code without due intervention of
the author<br>
</li>
<li>usually this object will be generated by a generating
function<br>
</li>
<li>the author will provide a couple of method-classes and
generating functions </li>
</ul>
</ul>
</li>
</ul>
</li>
</ul>
<h3>output structure</h3>
<ul>
<li> time-stamped elements: <br>
</li>
<ul>
<li>at least of type&nbsp;
<span style="font-family: Courier New,Courier,monospace;">ts</span>,
<br>
</li>
<li>if irregular input
<span style="font-family: Courier New,Courier,monospace;">zoo</span>
</li>
</ul>
<li> output type <br>
</li>
<ul>
<li>at least: S3-classes, <br>
</li>
</ul>
<ul>
<li>better:
S4-classes </li>
</ul>
<li> inheritance</li>
<ul>
<li>class should ``inherit'' from corresponding
classical return class <br>
</li>
<li>if return object of robust procedure is S4 and of the classical
one is S3:
<br>
S4-class should have the same slots as classical S3-class </li>
</ul>
<li> in case of several methods computed in parallel: <br>
</li>
<ul>
<li>list of corresponding
objects -&gt; <br>
</li>
<ul>
<li>particular plot/print methods for such
lists <br>
</li>
<li>particular consistency checks / validity methods<br>
</li>
</ul>
</ul>
<li> particular / overloaded methods
like <span style="font-family: Courier New,Courier,monospace;">print.acf</span>
or <span style="font-family: Courier New,Courier,monospace;">getMethod("print","acfrob")</span>
<ul>
<li> if possible, every argument legal for least specific
method <br>
(eg. <span style="font-family: Courier New,Courier,monospace;">plot.default</span>)
should be usable in particular (<span
style="font-family: Courier New,Courier,monospace;">plot</span>-)function
</li>
</ul>
</li>
</ul>

<br>
<hr style="width: 100%; height: 4px;">
<h1><br>
</h1>
<h2>Links to additional documents</h2>
<br>
<ul>
<li><big><a href="OpenIssues.txt">Open issues as to using R-Forge<br>
<br>
</a></big></li>
<li><big><a href="HOWTO-collaborate.txt">"HOWTO": What you have to do
to collaborate in 10 steps</a></big></li>
</ul>
<br>
<br>
<hr style="width: 100%; height: 4px;">

<h2>Maintained Target / Todo-List for Project "robust-ts"</h2>
<br>
<table style="text-align: left; width: 100%;" border="1" cellpadding="2"
cellspacing="2">
<tbody>
<tr>
<th style="vertical-align: top;">to be robustified: command from
package "stats"<br>
</th>
<th style="vertical-align: top;">robustification / reference<br>
</th>
<th style="vertical-align: top;">source of code<br>
</th>
<th style="vertical-align: top;">status<br>
</th>
<th style="vertical-align: top;">modus<br>
</th>
<th style="vertical-align: top;">working on it<br>
</th>
<th style="vertical-align: top;">from<br>
</th>
<th style="vertical-align: top;">to<br>
</th>
</tr>
<tr>
<td style="vertical-align: top;">acf / pacf<br>
</td>
<td style="vertical-align: top;">Ma/Genton<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">done upto user interface<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">now<br>
</td>
</tr>
<tr>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">quadrant correlation<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">M-estimator<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">AIC/BIC<br>
</td>
<td style="vertical-align: top;">Ronchetti<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">ar / arima<br>
</td>
<td style="vertical-align: top;">GM-estimators<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">tau-estimators, diagnostics<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">[g]arch </td>
<td style="vertical-align: top;">Boudt<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">code goes into fgarch package;<br>
here: just a wrapper<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">unitroot tests </td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">&nbsp;----<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">filter </td>
<td style="vertical-align: top;">see packages for Robust Kalman
Filtering / Robust Signal Extraction<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">Holt-winters </td>
<td style="vertical-align: top;">Gelper, Fried, Croux<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">Sarah Gelper<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">spec/spectrum </td>
<td style="vertical-align: top;">Spangl<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;">Bernhard Spangl<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">methods to be adapted:<br>
plot, print, summary / print.summary, confint, predict, residuals<br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;">to be added: function names from
Fin-Metrics </td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top;"><br>
</td>
</tr>
</tbody>
</table>
<br>
<h3>NOTE: </h3>This table is far from complete; so far it is just a starting point;
please feel free to complete the list; <br>
in particular fill in the exact references of papers; <br>
<br>
<h3>a sort of legend</h3>
<ul>
<li>column "source of code" is to capture the fact that maybe we get
code from outside, <br>
in particular from Insightful...</li>
<li>column "modus" will mention if original code resides in another R
package and is<br>
simply integrated by some wrapper --- as in the [g]arch case</li>
<li>more than one person may be listed in the "working on it"; the
"from" and "to"<br>
columns are meant to show in which time period a certain collaborator
has been/ was<br>
working on this topic<br>
</li>
</ul>
<hr style="width: 100%; height: 4px;">

<p> The <strong>project summary page</strong> you can find 
<a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>
here</strong></a>. </p>

</body>
</html>

