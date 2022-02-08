<html>
<head> 
<title></title> 
</head> 
<body> 
<h2></h2> 
 <?php 
 $path_to_folder="./"; 
 $dir=opendir($path_to_folder); 
 echo "<ul>";
 while($file=readdir($dir)) 
  { 
   if($file != "." && $file != ".." && $file != "index.php")
  {
   echo"<li><a href=\"$file\" target=\"blank\">$file</a></li>\n"; 
  }
 } 
 echo "</ul>";
 closedir($dir); 
?>
<?
##
?>
</body> 
</html>
