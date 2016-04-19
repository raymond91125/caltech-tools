#!/usr/bin/perl 

# flat file for genes is not on cronjob.  this machine probably does not have tazendra access, copied from .204  2016 03 24

use CGI;
use strict;
use LWP::Simple;
use LWP::UserAgent;
use Time::HiRes qw( time );

my $startTime = time; my $prevTime = time;
$startTime =~ s/(\....).*$/$1/;
$prevTime  =~ s/(\....).*$/$1/;

# use DBI;
# my $dbh = DBI->connect ( "dbi:Pg:dbname=testdb;host=131.215.52.76", "", "") or die "Cannot connect to database!\n";     # for remote access
# # my $dbh = DBI->connect ( "dbi:Pg:dbname=wobrdb", "", "") or die "Cannot connect to database!\n";
# my $result;

my $query = new CGI;
my $base_solr_url = 'http://wobr.caltech.edu:8082/solr/';		# raymond dev URL 2015 07 24

my $title = 'Tissue Enrichment Analysis';
my ($header, $footer) = &cshlNew($title);
&process();

sub process {
  my $action;                   # what user clicked
  unless ($action = $query->param('action')) { $action = 'none'; }

#   if ($action eq 'Tree') { &dag(); }
#     elsif ($action eq 'showGenes')                  { &showGenes();             }
#     elsif ($action eq 'queryChildren')              { &queryChildren();         }
#     elsif ($action eq 'annotSummaryGraph')          { &annotSummaryGraph();     }
#     elsif ($action eq 'annotSummaryCytoscape')      { &annotSummaryCytoscape(); }
#     elsif ($action eq 'annotSummaryJson')           { &annotSummaryJson();      }	# temporarily keep this for the live www.wormbase going through the fake phenotype_graph_json widget
#     elsif ($action eq 'annotSummaryJsonp')          { &annotSummaryJsonp();     }	# new jsonp widget to get directly from .wormbase without fake widget
  if ($action eq 'anatomySobaInput')           { &anatomySobaInput();      }
    elsif ($action eq 'HGE Analyze List')           { &anatomySoba('textarea'); }
    elsif ($action eq 'HGE Analyze File')           { &anatomySoba('file');     }
#     elsif ($action eq 'getWithPost')                { &getWithPost();           }	# test getting data with post, can delete later

#     elsif ($action eq 'showInferredGenes') { &showInferredGenes(); }	# combined into the single &showGenes();
#     elsif ($action eq 'showDirectGenes') {   &showDirectGenes();   }	# combined into the single &showGenes();
    else { &anatomySobaInput(); }				# no action, show dag by default
} # sub process

sub anatomySobaInput {
  &printHtmlHeader(); 
  print qq(<form method="post" action="tissue_enrichment_analysis.cgi" enctype="multipart/form-data">);
  print qq(Enter list of genes here<br/>);
  print qq(<textarea name="genelist" rows="20" cols="60"></textarea><br/>);
#   print qq(<input Type="checkbox" name="showProcessTimes" Value="showProcessTimes">Show Process Times<br/>\n);
  print qq(<input Type="checkbox" name="convertGeneToId" Value="convertGeneToId">Convert Genes to IDs<br/>\n);
#   print qq(<input Type="checkbox" name="calculateLcaNodes" Value="calculateLcaNodes">Calculate LCA Nodes<br/>\n);
  print qq(<input type="submit" name="action" value="HGE Analyze List"><br/><br/>\n);
  print qq(Upload a file with gene names :<br/>);
  print qq(<input type="file" name="geneNamesFile" /><br/>);
  print qq(<input type="submit" name="action" value="HGE Analyze File"><br/>\n);
  print qq(</form>);
#   &printMessageFooter(); 		# raymond wanted to remove this 2016 04 14
  &printHtmlFooter(); 
} # sub anatomySobaInput

sub anatomySoba {
  my ($filesource) = @_;
  &printHtmlHeader(); 

  my ($var, $datatype)          = &getHtmlVar($query, 'datatype');
#   ($var, my $showProcessTimes)  = &getHtmlVar($query, 'showProcessTimes');
  ($var, my $convertGeneToId)   = &getHtmlVar($query, 'convertGeneToId');
  ($var, my $calculateLcaNodes) = &getHtmlVar($query, 'calculateLcaNodes');
#   my ($var, $download)    = &getHtmlVar($query, 'download');
  unless ($datatype) { $datatype = 'anatomy'; }			# later will need to change based on different datatypes

#   if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Loading dictionary"); print qq($message<br/>\n); }

# python3  /home/raymond/local/src/git/TissueEnrichmentAnalysis/hypergeometricTests.py  /home/raymond/local/src/git/dictionary_generator/anat_dict.csv  /home/raymond/local/src/git/TissueEnrichmentAnalysis/input/SCohen_daf22.csv -p -s -t "BLAH"


  my %dict;
  my $dictFile = '/home/raymond/local/src/git/dictionary_generator/anat_dict.csv';
  open (DICT, "<$dictFile") or die "Cannot open $dictFile : $!";
  while (my $line = <DICT>) {
    my (@stuff) = split/,/, $line;
    if ($stuff[0] =~ m/WBGene/) { $dict{$stuff[0]}++; }
  } # while (my $line = <DICT>)
  close (DICT) or die "Cannot close $dictFile : $!";

#   if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Getting altId mappings"); print qq($message<br/>\n); }

  my @annotAnatomyTerms;					# array of annotated terms to loop and do pairwise comparisons
  my $genelist = '';
  if ($filesource eq 'textarea') {
      ($var, $genelist) = &getHtmlVar($query, 'genelist'); }
    elsif ($filesource eq 'file') {
      my $upload_filehandle = $query->upload("geneNamesFile");
      while ( <$upload_filehandle> ) { $genelist .= $_; }
    }
  if ($genelist =~ m/,/) { $genelist =~ s/,/ /g; }

#   if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Getting postgres gene name mappings"); print qq($message<br/>\n); }
  my %geneNameToId; my %geneIdToName;
  if ($convertGeneToId) {
#     my ($geneNameToIdHashref, $geneIdToNameHashref) = &populateGeneNamesFromPostgres();
    my ($geneNameToIdHashref, $geneIdToNameHashref) = &populateGeneNamesFromFlatfile();
    %geneNameToId        = %$geneNameToIdHashref;
    %geneIdToName        = %$geneIdToNameHashref; }

#   my %geneAnatomy; my %anatomyGene;

#   if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Processing user genes for validity"); print qq($message<br/>\n); }
  my (@names) = split/\s+/, $genelist;
  unless ($convertGeneToId) { foreach my $name (@names) { $geneNameToId{lc($name)} = $name; $geneIdToName{$name} = $name; } }
  my @invalidGene; my @nodataGene; my @goodGene; 
  foreach my $name (@names) {
    my ($lcname) = lc($name);
    my $wbgene = '';
    if ($geneNameToId{$lcname}) {
        $wbgene = $geneNameToId{$lcname}; 
        if ($dict{$wbgene}) { push @goodGene, $wbgene; }
          else { push @nodataGene, $wbgene; } }
      else { push @invalidGene, $name; }
  } # foreach my $name (@names)

  my %anatomyTerms;
  if (scalar @goodGene > 0) {
#       if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Processing hgf"); print qq($message<br/>\n); }
      my $time = time;
      my $tempfile     = '/tmp/hyperGeo/hyperGeo' . $time;
      my $tempOutFile  = '/tmp/hyperGeo/hyperGeo' . $time . '.txt';
      my $tempOutUrl   = '../data/hyperGeo/hyperGeo' . $time . '.txt';
      open (OUT, ">$tempOutFile") or die "Cannot open $tempOutFile : $!";
      my $tempImageUrl = '../data/hyperGeo/hyperGeo' . $time . '.sgv';
      open (TMP, ">$tempfile") or die "Cannot open $tempfile : $!";
#       print TMP qq(gene,reads\n);
      foreach my $gene (@goodGene) { print TMP qq($gene\n); }
      close (TMP) or die "Cannot close $tempfile : $!";
#       my $hyperData = `python /home/raymond/local/src/git/tissue_enrichment_tool_hypergeometric_test/src/hypergeometricTests.py /home/azurebrd/local/src/git/tissue_enrichment_tool_hypergeometric_test/genesets/WBPaper00013489_Ray_Enriched_WBbt:0006941_25`;
#       my $hyperData = `python /home/azurebrd/public_html/cgi-bin/hypergeometricTests.py /home/azurebrd/local/src/git/tissue_enrichment_tool_hypergeometric_test/genesets/WBPaper00013489_Ray_Enriched_WBbt:0006941_25`;
#       my $hyperData = `python /home/azurebrd/public_html/cgi-bin/hypergeometricTests.py $tempfile`;
#       my $hyperData = `python /home/raymond/local/src/git/tissue_enrichment_tool_hypergeometric_test/src/hypergeometricTests.py $tempfile`;
#                      python3  /home/raymond/local/src/git/TissueEnrichmentAnalysis/hypergeometricTests.py  /home/raymond/local/src/git/dictionary_generator/anat_dict.csv  /home/raymond/local/src/git/TissueEnrichmentAnalysis/input/SCohen_daf22.csv -p -s -t "BLAH"
      my $hyperData = `/home/raymond/local/src/git/TissueEnrichmentAnalysis/bin/tea  -d /home/raymond/local/src/git/dictionary_generator/anat_dict.csv  $tempfile "$tempfile" -p -s`;

#       `rm $tempfile`;
# print qq(HPD $hyperData HPD<br/>);
#       ($hyperData) = $hyperData =~ m/------------------------\n(.*?)------------------------/ms;
#       if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Processing hgf results to display"); print qq($message<br/>\n); }
      my (@hyperData) = split/\n/, $hyperData;
#       my $header = shift @hyperData;
#       my (@header) = split/,/, $header;
#       my $th = join"</th><th>", @header;
      if (scalar @hyperData > 0) {
#           print qq(<table><tr><th>$th</th></tr>\n);
        print qq(<table border="1" style="border-spacing: 0;">);  
        foreach my $line (@hyperData) {
          next if ($line =~ m/Executing script/);
          if ($line) {
            print OUT qq($line\n);
            $line =~ s|\t|</td><td>|g;
            if ($line =~ m/(WBbt:\d+)/) { 
              my $wbbt = $1;
              my $url = 'http://www.wormbase.org/species/all/anatomy_term/' .$wbbt . '#013--10';
              $line =~ s/$wbbt/<a href="$url" target="_blank">$wbbt<\/a>/; }
            print qq(<tr><td align="right">$line</td></tr>);
          }
        } # foreach my $line (@hyperData)
        print qq(</table>);  

#           my %sort;
#           foreach my $line (@hyperData) {
#             my ($name, $wbbt, $value) = $line =~ m/^(.*?)\((WBbt:.*?)\),([^,]*?)$/;
#             my $url = 'http://www.wormbase.org/species/all/anatomy_term/' .$wbbt . '#013--10';
#             $sort{$value}{$name} = $url;
#           } # foreach my $line (@hyperData)
#           foreach my $value (sort {$a<=>$b} keys %sort) {
#             foreach my $name (sort keys %{ $sort{$value} }) {
#               my $url  = $sort{$value}{$name};
#               my $link = qq(<a href="$url" target="_blank">$name</a>);
# #               print qq($name,$value<br/>);
#               print qq(<tr><td>$link</td><td>$value</td></tr>);
#             } # foreach my $name (sort keys %{ $sort{$value} })
#           } # foreach my $value (sort keys %sort)
#           print qq(</table>\n);
        }
        else { print qq(No significantly enriched cell/tissue has been found.<br/>\n); }
      print qq(<img src="$tempImageUrl"><br/>Drag graph to your desktop to save.<br/>);
      close (OUT) or die "Cannot close $tempOutFile : $!";
      print qq(Download output table <a href="$tempOutUrl" target="_blank">here</a><br/><br/>);
    }
    else { print qq(There are no genes with annotated data to generate results.<br/>\n); }
  print qq(<br/>);
#   if ($showProcessTimes) { (my $message) = &getDiffTime($startTime, $prevTime, "Displaying gene sets"); print qq($message<br/>\n); }
  if (scalar @invalidGene > 0) {
    my $countInvalidGenes = scalar @invalidGene;
    print qq(Your list has $countInvalidGenes invalid WormBase genes :<br/>\n);
    print qq(<textarea rows="6" cols="80">);
    foreach my $gene (@invalidGene) { print qq($gene\n); } 
    print qq(</textarea><br/><br/>); }
  if (scalar @nodataGene > 0) {
    my $countNodataGenes = scalar @nodataGene;
    print qq(Your list has $countNodataGenes valid WormBase genes that have no annotated data or are excluded from testing :<br/>\n);
    print qq(<textarea rows="6" cols="80">);
    foreach my $gene (@nodataGene) { print qq($gene - $geneIdToName{$gene}\n); } 
    print qq(</textarea><br/><br/>); }
  if (scalar @goodGene > 0) {
    my $countGoodGenes = scalar @goodGene;
    print qq(Your list has $countGoodGenes valid WormBase genes included in statistical testing :<br/>\n);
    print qq(<textarea rows="6" cols="80">);
    foreach my $gene (@goodGene) { print qq($gene - $geneIdToName{$gene}\n); } 
    print qq(</textarea><br/><br/>); }
  print qq(<a href="tissue_enrichment_analysis.cgi">perform another query</a><br/>);

#   &printMessageFooter(); 		# raymond wanted to remove this 2016 04 14
  &printHtmlFooter(); 
# http://131.215.12.204/~azurebrd/cgi-bin/amigo.cgi?genelist=WBGene00010209+WBGene00010212+WBGene00010295+WBGene00015814+WBGene11111111+WBGene00000012+WBGene00000013+WBGene00000002+ZK512.6%0D%0A&action=anatomySoba
# http://131.215.12.204/~azurebrd/cgi-bin/amigo.cgi?action=anatomySobaInput
# WBGene00010209 WBGene00010212 WBGene00010295 WBGene00015814 WBGene11111111 WBGene00000012 WBGene00000013 WBGene00000002 ZK512.6 qwera 2z3zt daf-2
} # sub anatomySoba


sub populateGeneNamesFromFlatfile {
  my %geneNameToId; my %geneIdToName;
  my $infile = '/home/azurebrd/cron/gin_names/gin_names.txt';
  open (IN, "<$infile") or die "Cannot open $infile : $!";
  while (my $line = <IN>) {
    chomp $line;
    my ($id, $name, $primary) = split/\t/, $line;
    if ($primary eq 'primary') { $geneIdToName{$id}     = $name; }
    my ($lcname)           = lc($name);
    $geneNameToId{$lcname} = $id; }
  close (IN) or die "Cannot close $infile : $!";
  return (\%geneNameToId, \%geneIdToName);
} # sub populateGeneNamesFromFlatfile


sub printMessageFooter { print qq(if you use this tool, please cite us. If you have ideas, requests or bug-issues, please contact <a href="mailto:dangeles\@caltech.edu">dangeles\@caltech.edu</a><br/>); }

# sub printHtmlFooter { print qq(</body></html>\n); }
sub printHtmlFooter { print $footer }

sub printHtmlHeader { 
  my $javascript = << "EndOfText";
<script src="http://code.jquery.com/jquery-1.9.1.js"></script>
<script type="text/javascript">
function toggleShowHide(element) {
    document.getElementById(element).style.display = (document.getElementById(element).style.display == "none") ? "" : "none";
    return false;
}
function togglePlusMinus(element) {
    document.getElementById(element).innerHTML = (document.getElementById(element).innerHTML == "&nbsp;+&nbsp;") ? "&nbsp;-&nbsp;" : "&nbsp;+&nbsp;";
    return false;
}
</script>
EndOfText
#   print qq(Content-type: text/html\n\n<html><head><title>Tissue Enrichment Analysis</title>$javascript</head><body>\n); 
  print qq(Content-type: text/html\n\n$header $javascript<body>\n); 
}

sub getHtmlVar {                
  no strict 'refs';             
  my ($query, $var, $err) = @_; 
  unless ($query->param("$var")) {
    if ($err) { print "<FONT COLOR=blue>ERROR : No such variable : $var</FONT><BR>\n"; }
  } else { 
    my $oop = $query->param("$var");
    $$var = &untaint($oop);         
    return ($var, $$var);           
  } 
} # sub getHtmlVar

sub untaint {
  my $tainted = shift;
  my $untainted;
  if ($tainted eq "") {
    $untainted = "";
  } else { # if ($tainted eq "")
    $tainted =~ s/[^\w\-.,;:?\/\\@#\$\%\^&*\>\<(){}[\]+=!~|' \t\n\r\f\"€‚ƒ„…†‡ˆ‰Š‹ŒŽ‘’“”•—˜™š›œžŸ¡¢£¤¥¦§¨©ª«¬­®¯°±²³´µ¶·¹º»¼½¾¿ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖ×ØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõö÷øùúûüýþ]//g;
    if ($tainted =~ m/^([\w\-.,;:?\/\\@#\$\%&\^*\>\<(){}[\]+=!~|' \t\n\r\f\"€‚ƒ„…†‡ˆ‰Š‹ŒŽ‘’“”•—˜™š›œžŸ¡¢£¤¥¦§¨©ª«¬­®¯°±²³´µ¶·¹º»¼½¾¿ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖ×ØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõö÷øùúûüýþ]+)$/) {
      $untainted = $1;
    } else {
      die "Bad data Tainted in $tainted";
    }
  } # else # if ($tainted eq "")
  return $untainted;
} # sub untaint

sub cshlNew {
  my $title = shift;
  unless ($title) { $title = ''; }      # init title in case blank
  my $page = get "http://tazendra.caltech.edu/~azurebrd/sanger/wormbaseheader/WB_header_footer.html";
#  $page =~ s/href="\//href="http:\/\/www.wormbase.org\//g;
#  $page =~ s/src="/src="http:\/\/www.wormbase.org/g;
  ($header, $footer) = $page =~ m/^(.*?)\s+DIVIDER\s+(.*?)$/s;  # 2006 11 20    # get this from tazendra's script result.
#   $header =~ s/WormBase - Home Page/$title/g;                 # 2015 05 07    # wormbase 2.0
  $header =~ s/<title>.*?<\/title>/<title>$title<\/title>/g;
  return ($header, $footer);
} # sub cshlNew

