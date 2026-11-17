param(
  [ValidateRange(1,8)][int]$Workers=2,
  [switch]$KeepAlive,
  [string]$Resume='',
  [string]$RscriptPath='E:/应用/R-4.6.1/bin/Rscript.exe'
)
$ErrorActionPreference='Stop'
$project=[IO.Directory]::GetParent($PSScriptRoot).Parent.FullName
$entry=Join-Path $PSScriptRoot 'joint_TS7_sensitivity.R'
if(-not(Test-Path -LiteralPath $RscriptPath)){throw "Rscript not found: $RscriptPath"}
$env:LC_ALL='English_United States.utf8'
$env:LANG='en_US.UTF-8'
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$checkArgs=@('--vanilla','--encoding=UTF-8',$entry,"--project=$project","--workers=$Workers")
if($Resume){$checkArgs+="--resume=$Resume"}
& $RscriptPath @checkArgs
if($LASTEXITCODE -ne 0){throw 'Preflight failed; background job was not started.'}
$logParent=Join-Path $project 'archive/joint_TS7_sensitivity_launch'
New-Item -ItemType Directory -Path $logParent -Force | Out-Null
$logs=$null
for($i=1;$i -le 9999;$i++){
  $candidate=Join-Path $logParent ('launch{0:d3}' -f $i)
  if(-not(Test-Path -LiteralPath $candidate)){$logs=(New-Item -ItemType Directory -Path $candidate).FullName;break}
}
if(-not $logs){throw 'Cannot create launch log directory.'}
$arguments=@('--vanilla','--encoding=UTF-8',('"'+$entry+'"'),'--run',('"--project='+$project+'"'),"--workers=$Workers")
if($Resume){$arguments+=('"--resume='+$Resume+'"')}
$child=Start-Process -FilePath $RscriptPath -ArgumentList $arguments -WorkingDirectory $project -WindowStyle Hidden -PassThru `
  -RedirectStandardOutput (Join-Path $logs 'stdout.log') -RedirectStandardError (Join-Path $logs 'stderr.log')
$record=[ordered]@{pid=$child.Id;workers=$Workers;project=$project;resume=$Resume;logs=$logs;started=(Get-Date).ToString('o')}
$record | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $logs 'launch.json') -Encoding UTF8
$record | ConvertTo-Json
if($KeepAlive){Wait-Process -Id $child.Id -ErrorAction SilentlyContinue}
