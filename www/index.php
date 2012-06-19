
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

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
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

  <hr>
  <b>Description:</b>
  
  The R-package <tt>hypergsplines</tt> implements objective Bayesian variable
  selection in generalized additive models (Normal, Binomial and Poisson). As
  objective parameters prior a (generalized) hyper-g prior is used. As objective
  model prior so-called multiplicity correction priors are used. More details on
  the methodology can be found in <a href="http://arxiv.org/abs/1108.3520"> the
    technical report (2011)</a>.
  <hr>
  <b>Download:</b>

  Currently the R-package <tt>hypergsplines</tt> is not available from CRAN,
  because it is still under continuing development. The current versions are
  available from R-Forge (this site). If you are running a recent R version, you
  can obtain a binary snapshot of the development version as 
  <div style="text-align:
	      center;"> 
    <tt>install.packages("hypergsplines",repos="http://r-forge.r-project.org")</tt>.
  </div>
  You can also manually download
  the <a href="http://r-forge.r-project.org/bin">binary snapshot</a> or
  the <a href="http://r-forge.r-projec.org/src">source tarball</a>. 
  You can have a look at the recent changes to the
  package <a href="https://r-forge.r-project.org/scm/viewvc.php/pkg/ChangeLog?view=markup&root=hypergsplines">here</a>.   
  The project summary page can be found <a href="http://<?php echo $domain;
  ?>/projects/<?php echo $group_name; ?>/">here</a>. 
  <hr>
  <b>Documentation:</b>

  The R-package features a vignette which introduces the most important functions
  in a logistic additive regression example. After installing the package, you can
  open the vignette from within R with the command
  <div style="text-align:
	      center;"> 
    <tt>vignette("examples", package="hypergsplines")</tt>.
  </div>
  You can also have a look at the ISBA 2012 <a href="poster.pdf">poster</a>
  presenting the methodology. 
  <hr>
  <b>Developers:</b>

  <a href="http://www.biostat.uzh.ch/aboutus/people/sabanes.html">Daniel Sabanes
  Bove</a>, Division of Biostatistics, Institute of Social and Preventive
  Medicine, University of Zurich, Switzerland

</body>
</html>
