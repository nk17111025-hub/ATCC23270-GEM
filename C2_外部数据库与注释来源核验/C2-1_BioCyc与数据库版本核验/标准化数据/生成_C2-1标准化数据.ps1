$root = 'D:\嗜酸氧化亚铁硫杆菌\03_BioCyc\C2-1'
$raw = Join-Path $root '官方原始下载'
$std = Join-Path $root '标准化数据'

function C($v) {
  if ($null -eq $v) { return '' }
  return (([string]$v) -replace '[\t\r\n]+',' ' -replace '\s{2,}',' ').Trim()
}
function TextOf($n,$xp) {
  $q = $n.SelectSingleNode($xp)
  if ($null -eq $q) { return '' }
  return C $q.InnerText
}
function Vals($n,$xp) {
  $a = @($n.SelectNodes($xp) | ForEach-Object { C $_.InnerText } | Where-Object { $_ } | Select-Object -Unique)
  return ($a -join ';')
}
function Refs($n,$xp) {
  $a = @($n.SelectNodes($xp) | ForEach-Object {
    if ($_.Attributes['frameid']) { C $_.Attributes['frameid'].Value }
    elseif ($_.Attributes['ID']) { C $_.Attributes['ID'].Value }
  } | Where-Object { $_ } | Select-Object -Unique)
  return ($a -join ';')
}
function Links($n) {
  $a = @($n.SelectNodes('.//dblink') | ForEach-Object {
    $db = TextOf $_ './dblink-db'; $id = TextOf $_ './dblink-oid'
    if ($db -and $id) { '{0}:{1}' -f $db,$id }
  } | Where-Object { $_ } | Select-Object -Unique)
  return ($a -join ';')
}
function LoadLatest($pat) {
  $p = Get-ChildItem -LiteralPath $raw -Filter $pat | Sort-Object LastWriteTime | Select-Object -Last 1
  $x = New-Object System.Xml.XmlDocument
  $x.Load($p.FullName)
  return [pscustomobject]@{ Path=$p.FullName; Xml=$x }
}
function WriteTsv($path,$headers,$rows) {
  $lines = @($headers -join "`t")
  foreach ($r in $rows) {
    $v = @()
    foreach ($h in $headers) { $v += C $r[$h] }
    $lines += ($v -join "`t")
  }
  [IO.File]::WriteAllLines($path,$lines,[Text.UTF8Encoding]::new($false))
}

$g = LoadLatest '*all-genes_full_*.xml'
$p = LoadLatest '*proteins_full_retry_*.xml'
$rx = LoadLatest 'GCF_000021485_reactions_full_*.xml'
$pw = LoadLatest '*pathways_full_*.xml'
$er = LoadLatest '*enzymatic-reactions_full_*.xml'
$compoundFiles = @(Get-ChildItem -LiteralPath $raw -Filter '*compounds_*network_full_*.xml' | Sort-Object LastWriteTime)

$genes=@(); $proteins=@(); $reactions=@(); $pathways=@(); $enzymes=@(); $relations=@(); $evidence=@(); $compounds=@{}

foreach ($n in @($g.Xml.DocumentElement.SelectNodes('./Gene'))) {
  $id=$n.GetAttribute('frameid'); $prod=Refs $n './product/*'
  $genes += [ordered]@{类别='Gene';PGDB='GCF_000021485';对象ID=$id;名称=(Vals $n './common-name|./synonym');同义词=(Vals $n './synonym');accession=(Vals $n './accession-1|./accession-2');product=$prod;parent=(Refs $n './parent/*');replicon=(Refs $n './replicon/*');方向=(Vals $n './transcription-direction');左端=(Vals $n './left-end-position');右端=(Vals $n './right-end-position');来源=$g.Path}
  foreach ($q in @($n.SelectNodes('./product/*'))) { $relations += [ordered]@{关系类型='gene_to_protein';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$g.Path;备注='gene/product'} }
}

foreach ($n in @($p.Xml.DocumentElement.SelectNodes('./Protein'))) {
  $id=$n.GetAttribute('frameid'); $genesRef=Refs $n './gene/*'
  $proteins += [ordered]@{类别='Protein';PGDB='GCF_000021485';对象ID=$id;名称=(Vals $n './common-name');基因=$genesRef;注释=(Vals $n './comment');分子量_kD=(Vals $n './molecular-weight-seq');外部ID=(Links $n);来源=$p.Path}
  foreach ($q in @($n.SelectNodes('./gene/*'))) { $relations += [ordered]@{关系类型='protein_to_gene';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$p.Path;备注='protein/gene'} }
}

foreach ($n in @($rx.Xml.DocumentElement.SelectNodes('./Reaction'))) {
  $id=$n.GetAttribute('frameid'); $left=Refs $n './left/*'; $right=Refs $n './right/*'; $enz=Refs $n './/enzyme/*'; $enzrx=Refs $n './enzymatic-reaction/*'; $pwys=Refs $n './in-pathway/*'
  $reactions += [ordered]@{类别='Reaction';PGDB='GCF_000021485';对象ID=$id;名称=(Vals $n './common-name');方向=(Vals $n './reaction-direction');左侧=$left;右侧=$right;酶=$enz;酶促反应=$enzrx;通路=$pwys;EC号=(Vals $n './ec-number');孤儿=(Vals $n './orphan');生理相关=(Vals $n './physiologically-relevant');证据代码=(Refs $n './/Evidence-Code');文献=(Refs $n './/Publication');外部ID=(Links $n);来源=$rx.Path}
  foreach ($q in @($n.SelectNodes('./left/*'))) { $relations += [ordered]@{关系类型='reaction_left';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$rx.Path;备注=$q.Name} }
  foreach ($q in @($n.SelectNodes('./right/*'))) { $relations += [ordered]@{关系类型='reaction_right';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$rx.Path;备注=$q.Name} }
  foreach ($q in @($n.SelectNodes('.//enzyme/*'))) { $relations += [ordered]@{关系类型='reaction_to_enzyme';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$rx.Path;备注='enzymatic-reaction/enzyme'} }
  foreach ($q in @($n.SelectNodes('./enzymatic-reaction/*'))) { $relations += [ordered]@{关系类型='reaction_to_enzymatic_reaction';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$rx.Path;备注='reaction/enzymatic-reaction'} }
  foreach ($q in @($n.SelectNodes('./in-pathway/*'))) { $relations += [ordered]@{关系类型='reaction_to_pathway';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$rx.Path;备注='in-pathway'} }
  foreach ($q in @($n.SelectNodes('.//Evidence-Code'))) { $evidence += [ordered]@{主体类型='Reaction';主体ID=$id;证据类型='Evidence-Code';证据ID=$q.GetAttribute('frameid');名称=(TextOf $q './common-name');文献ID='';来源=$rx.Path} }
  foreach ($q in @($n.SelectNodes('.//Publication'))) { $evidence += [ordered]@{主体类型='Reaction';主体ID=$id;证据类型='Publication';证据ID=$q.GetAttribute('frameid');名称=(TextOf $q './title');文献ID=$q.GetAttribute('frameid');来源=$rx.Path} }
}

foreach ($n in @($pw.Xml.DocumentElement.SelectNodes('./Pathway'))) {
  $id=$n.GetAttribute('frameid')
  $pathways += [ordered]@{类别='Pathway';PGDB='GCF_000021485';对象ID=$id;名称=(Vals $n './common-name');反应=(Refs $n './/Reaction');子通路=(Refs $n './sub-pathways/*');基因=(Refs $n './/Gene');证据代码=(Refs $n './/Evidence-Code');文献=(Refs $n './/Publication');来源=$pw.Path}
  foreach ($q in @($n.SelectNodes('.//Reaction'))) { $relations += [ordered]@{关系类型='pathway_to_reaction';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$pw.Path;备注='pathway reaction'} }
  foreach ($q in @($n.SelectNodes('./sub-pathways/*'))) { $relations += [ordered]@{关系类型='pathway_to_subpathway';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$pw.Path;备注='sub-pathway'} }
  foreach ($q in @($n.SelectNodes('.//Gene'))) { $relations += [ordered]@{关系类型='pathway_to_gene';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$pw.Path;备注='pathway gene'} }
  foreach ($q in @($n.SelectNodes('.//Evidence-Code'))) { $evidence += [ordered]@{主体类型='Pathway';主体ID=$id;证据类型='Evidence-Code';证据ID=$q.GetAttribute('frameid');名称=(TextOf $q './common-name');文献ID='';来源=$pw.Path} }
  foreach ($q in @($n.SelectNodes('.//Publication'))) { $evidence += [ordered]@{主体类型='Pathway';主体ID=$id;证据类型='Publication';证据ID=$q.GetAttribute('frameid');名称=(TextOf $q './title');文献ID=$q.GetAttribute('frameid');来源=$pw.Path} }
}

foreach ($n in @($er.Xml.DocumentElement.SelectNodes('./Enzymatic-Reaction'))) {
  $id=$n.GetAttribute('frameid')
  $enzymes += [ordered]@{类别='Enzymatic-Reaction';PGDB='GCF_000021485';对象ID=$id;名称=(Vals $n './common-name');反应=(Refs $n './reaction/*');酶=(Refs $n './enzyme/*');EC号=(Vals $n './ec-number');来源=$er.Path}
  foreach ($q in @($n.SelectNodes('./reaction/*'))) { $relations += [ordered]@{关系类型='enzymatic_reaction_to_reaction';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$er.Path;备注='enzymatic reaction'} }
  foreach ($q in @($n.SelectNodes('./enzyme/*'))) { $relations += [ordered]@{关系类型='enzymatic_reaction_to_enzyme';主体ID=$id;客体ID=$q.GetAttribute('frameid');来源=$er.Path;备注='enzyme'} }
}

foreach ($cf in $compoundFiles) {
  $x=New-Object System.Xml.XmlDocument; $x.Load($cf.FullName)
  foreach ($n in @($x.DocumentElement.SelectNodes('./Compound'))) {
    $id=$n.GetAttribute('frameid'); if (-not $id) { continue }
    if (-not $compounds.ContainsKey($id)) { $compounds[$id]=[ordered]@{类别='Compound';PGDB='GCF_000021485';对象ID=$id;名称=(Vals $n './common-name');同义词=(Vals $n './synonym');化学式=(Vals $n './chemical-formula');分子量=(Vals $n './molecular-weight');外部ID=(Links $n);来源=$cf.FullName} }
  }
}

WriteTsv (Join-Path $std 'genes.tsv') @('类别','PGDB','对象ID','名称','同义词','accession','product','parent','replicon','方向','左端','右端','来源') $genes
WriteTsv (Join-Path $std 'proteins.tsv') @('类别','PGDB','对象ID','名称','基因','注释','分子量_kD','外部ID','来源') $proteins
WriteTsv (Join-Path $std 'reactions.tsv') @('类别','PGDB','对象ID','名称','方向','左侧','右侧','酶','酶促反应','通路','EC号','孤儿','生理相关','证据代码','文献','外部ID','来源') $reactions
WriteTsv (Join-Path $std 'pathways.tsv') @('类别','PGDB','对象ID','名称','反应','子通路','基因','证据代码','文献','来源') $pathways
WriteTsv (Join-Path $std 'enzymatic-reactions.tsv') @('类别','PGDB','对象ID','名称','反应','酶','EC号','来源') $enzymes
WriteTsv (Join-Path $std 'compounds_network.tsv') @('类别','PGDB','对象ID','名称','同义词','化学式','分子量','外部ID','来源') @($compounds.Values)
WriteTsv (Join-Path $std 'relations.tsv') @('关系类型','主体ID','客体ID','来源','备注') $relations
WriteTsv (Join-Path $std 'evidence_references.tsv') @('主体类型','主体ID','证据类型','证据ID','名称','文献ID','来源') $evidence
WriteTsv (Join-Path $root 'C2-1_基因_蛋白_反应_通路.tsv') @('关系类型','主体ID','客体ID','来源','备注') $relations
WriteTsv (Join-Path $root 'C2-1_文献与证据关联.tsv') @('主体类型','主体ID','证据类型','证据ID','名称','文献ID','来源') $evidence
Write-Output ("genes={0} proteins={1} reactions={2} pathways={3} enzymatic_reactions={4} compounds_unique={5} relations={6} evidence={7}" -f $genes.Count,$proteins.Count,$reactions.Count,$pathways.Count,$enzymes.Count,$compounds.Count,$relations.Count,$evidence.Count)
