param([ValidateRange(100,100000)][int]$Draws=4000,
      [string]$RscriptPath='E:/应用/R-4.6.1/bin/Rscript.exe')
$ErrorActionPreference='Stop'
$project=[IO.Directory]::GetParent($PSScriptRoot).Parent.FullName
$entry=Join-Path $PSScriptRoot 'joint_TS7_prior_predictive.R'
$env:LC_ALL='English_United States.utf8'
$env:LANG='en_US.UTF-8'
$env:OMP_NUM_THREADS='1'
& $RscriptPath --vanilla --encoding=UTF-8 $entry --self-test
if($LASTEXITCODE -ne 0){throw 'Distribution self-test failed'}
& $RscriptPath --vanilla --encoding=UTF-8 $entry "--project=$project" "--draws=$Draws"
if($LASTEXITCODE -ne 0){throw 'Input check failed'}
$parent=Join-Path $project 'archive/joint_TS7_prior_predictive_launch'
New-Item -ItemType Directory -Path $parent -Force | Out-Null
$logs=$null
for($i=1;$i -le 9999;$i++){
  $candidate=Join-Path $parent ('launch{0:d3}' -f $i)
  if(-not(Test-Path -LiteralPath $candidate)){$logs=(New-Item -ItemType Directory -Path $candidate).FullName;break}
}
if(-not $logs){throw 'No available log directory'}
$arguments=@('--vanilla','--encoding=UTF-8',('"'+$entry+'"'),'--run',('"--project='+$project+'"'),"--draws=$Draws")
$child=Start-Process -FilePath $RscriptPath -ArgumentList $arguments -WorkingDirectory $project -WindowStyle Hidden -PassThru `
  -RedirectStandardOutput (Join-Path $logs 'stdout.log') -RedirectStandardError (Join-Path $logs 'stderr.log')
$record=[ordered]@{pid=$child.Id;draws=$Draws;logs=$logs;project=$project;started=(Get-Date).ToString('o')}
$record | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $logs 'launch.json') -Encoding UTF8
$record | ConvertTo-Json
