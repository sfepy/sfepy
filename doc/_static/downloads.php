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
      <a href="https://sfepy.org"><img src="sfepy_logo_title_small.png" border="0" height="70px" alt="SfePy"/></a>
      <a href="https://ntc.zcu.cz/en/"><img src="ntc_logo.png" align="right" border="0" height="80px" alt="NTC"/></a>
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
  $con = mysqli_connect('127.0.0.1:3307', 'sfepy', 'sfepyheslo2018', 'sfepy');
  if (mysqli_connect_errno($con)) {
    die('Failed to connect to MySQL: ' . mysqli_connect_error());
  }

  echo '<table width="100%">';
  echo '<tr><td><i>Filename</i></td><td><i>UploadDate</i></td><td><i>Size</i></td><td><i>DownloadCount</i></td><td></td></tr>';
  // echo '<tr><td><i>Filename</i></td><td><i>UploadDate</i></td><td><i>Size</i></td><td></td></tr>';
  $tgz_files = glob($GLOBALS['download_dir'].'/*.{tar.gz}', GLOB_BRACE);
  $itgz_files = array_reverse($tgz_files, true);
  foreach($itgz_files as $file_) {
    $file = basename($file_);
    $sql = 'SELECT version,num_downloads FROM sfepy_downloads WHERE version="'.$file.'"';
    if (!($data = mysqli_query($con, $sql))) {
      die('DB Error, could not query the database: ' . mysqli_error($con));
    }

    if (mysqli_num_rows($data) == 0) {
      $sql = 'INSERT INTO sfepy_downloads values ("'.$file.'", 0)';
      mysqli_query($con, $sql);
      $version = $file;
      $num_downloads = 0;
    }
    else {
      list($version, $num_downloads) = mysqli_fetch_row($data);
    }

    $size = filesize($file_);
    $modtime = date("M d Y", filemtime($file_));
    $tag = '<tr><td><b>'.$version.'</b></td><td>'.$modtime.'</td><td>'.$size.' bytes</td><td>'.$num_downloads.'</td><td><a href="'.$_SERVER['PHP_SELF'].'?fun=download&ver='.$version.'"><button>Download</button></a></td></tr>';
    echo $tag;

    mysqli_free_result($data);
  }
  echo '</table>';

  mysqli_close($con);
}

function download($version) {
  $version_ = $GLOBALS['download_dir'].$version;
  header('Content-Description: File Transfer');
  header('Content-Type: application/octet-stream');
  header('Content-Disposition: attachment; filename='.basename($version_));
  header('Content-Transfer-Encoding: binary');
  header('Expires: 0');
  header('Cache-Control: must-revalidate');
  header('Pragma: public');
  header('Content-Length: ' . filesize($version_));
  ob_clean();
  flush();
  readfile($version_);

  $con = mysqli_connect('127.0.0.1:3307', 'sfepy', 'sfepyheslo2018', 'sfepy');
  if (mysqli_connect_errno($con)) {
    die('Failed to connect to MySQL: ' . mysqli_connect_error());
  }

  $sql = 'SELECT num_downloads FROM sfepy_downloads WHERE version="'.$version.'"';
  if (!($data = mysqli_query($con, $sql))) {
    die('DB Error, could not query the database: ' . mysqli_error($con));
  }

  list($num_downloads) = mysqli_fetch_row($data);
  $num_downloads = $num_downloads + 1;
  $sql = 'UPDATE sfepy_downloads SET num_downloads="'.$num_downloads.'" WHERE version="'.$version.'"';
  if (!mysqli_query($con, $sql)) {
    die('DB Error, could not query the database: ' . mysqli_error($con));
  }

  mysqli_free_result($data);
  mysqli_close($con);
}

function dwl($version) {
  $tag = $_SERVER['PHP_SELF'].'?fun=download2&ver='.$version;
  header('Refresh: 0; '.$tag);
  flush();
  $tag = 'https://sfepy.org/doc-devel/installation.html';
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
        &copy; Copyright 2018, Robert Cimrman and Contributors.
      Created using <a href="https://sphinx.pocoo.org/">Sphinx</a>.
    </div>
  </body>
</html>
