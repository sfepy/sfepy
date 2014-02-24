<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <title>SfePy: Simple Finite Elements in Python</title>

    <link rel="stylesheet" href="sfepy.css" type="text/css" />
    <link rel="stylesheet" href="pygments.css" type="text/css" />
  </head>
  <body>
    <div style="background-color: white; text-align: left; padding: 10px 10px 15px 15px">
      <a href="http://sfepy.org"><img src="sfepy_logo_title_small.png" border="0" height="70px" alt="SfePy"/></a>
      <a href="http://www.zcu.cz/ntc/en/"><img src="ntc_logo.png" align="right" border="0" height="80px" alt="NTC"/></a>
    </div>

    <div class="related">
      <h3>Navigation</h3>
      <ul>
      </ul>
    </div>

    <div class="document">
      <div class="body">
        <div class="section" id="sfepy-simple-finite-elements-in-python">
          <h1>SfePy: Downloads</h1>

          It is recommended to download the latest release of SfePy.<br><br>

<?php
$download_dir = '../../sfepy_downloads/';

function init(){
  $relace = mysql_connect('localhost', 'host', 'h0st207p');
  if (!$relace) {
    die('Could not connect: ' . mysql_error());
  }

  if (!mysql_select_db('host', $relace)) {
    die('Could not select database');
  }

  echo '<table width="100%">';
  echo '<tr><td><i>Filename</i></td><td><i>UploadDate</i></td><td><i>Size</i></td><td><i>DownloadCount</i></td><td></td></tr>';

  $tgz_files = glob($GLOBALS['download_dir'].'/*.{tar.gz}', GLOB_BRACE);
  $itgz_files = array_reverse($tgz_files, true);
  foreach($itgz_files as $file_) {
    $file = basename($file_);
    $sql = 'SELECT verze,stazeno FROM sfepy_downloads WHERE verze="'.$file.'"';
    $data = mysql_query($sql, $relace);

    if (!$data) {
      die('DB Error, could not query the database: ' . mysql_error());
    }

    if (mysql_num_rows($data) == 0) {
      $sql = 'INSERT INTO sfepy_downloads values ("'.$file.'", 0)';
      mysql_query($sql, $relace);
    }
    else {
      list($verze, $stazeno) = mysql_fetch_row($data);
      $size = filesize($file_);
      $modtime = date("M d Y", filemtime($file_));
      $tag = '<tr><td><b>'.$verze.'</b></td><td>'.$modtime.'</td><td>'.$size.' bytes</td><td>'.$stazeno.'</td><td><a href="'.$_SERVER['PHP_SELF'].'?fun=download&ver='.$verze.'"><button>Download</button></a></td></tr>';
      echo $tag;
    }

    mysql_free_result($data);
  }
  echo '</table>';

}

function download($verze) {
  $verze_ = $GLOBALS['download_dir'].$verze;
  header('Content-Description: File Transfer');
  header('Content-Type: application/octet-stream');
  header('Content-Disposition: attachment; filename='.basename($verze_));
  header('Content-Transfer-Encoding: binary');
  header('Expires: 0');
  header('Cache-Control: must-revalidate');
  header('Pragma: public');
  header('Content-Length: ' . filesize($verze_));
  ob_clean();
  flush();
  readfile($verze_);

  $relace = mysql_connect('localhost', 'host', 'h0st207p');
  if (!$relace) {
    die('Could not connect: ' . mysql_error());
  }

  if (!mysql_select_db('host', $relace)) {
    die('Could not select database');
  }

  $sql = 'SELECT stazeno FROM sfepy_downloads WHERE verze="'.$verze.'"';
  $data = mysql_query($sql, $relace);

  if (!$data) {
    die('DB Error, could not query the database: ' . mysql_error());
  }

  list($stazeno) = mysql_fetch_row($data);
  $stazeno = $stazeno + 1;
  $sql = 'UPDATE sfepy_downloads SET stazeno="'.$stazeno.'" WHERE verze="'.$verze.'"';
  mysql_query($sql, $relace);

  mysql_free_result($data);
}

function dwl($verze) {
  $tag = $_SERVER['PHP_SELF'].'?fun=download2&ver='.$verze;
  header('Refresh: 0; '.$tag);
  flush();
  $tag = 'http://sfepy.org/doc-devel/installation.html';
  echo '<br>Thank you for downloading SfePy! Check out the <a href="'.$tag.'">installation instructions.</a>';
}

if(!isset($_GET['fun'])){
  $fun = 'list';
} else {
  $fun = $_GET['fun'];
}

switch($fun) {
case 'list': init(); break;
case 'download2': download($_GET['ver']); break;
case 'download': dwl($_GET['ver']); break;
}
?>

        </div>
      </div>
    </div>
    <div class="related">
      <h3>Navigation</h3>
      <ul>
      </ul>
    </div>

    <div class="footer">
        &copy; Copyright 2010, Robert Cimrman and Contributors.
      Created using <a href="http://sphinx.pocoo.org/">Sphinx</a>.
    </div>
  </body>
</html>
