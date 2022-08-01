$source = 'http://rest.kegg.jp/get/br:ko00001/json'
$destination = './data/ko.json'

if (Test-Path $destination) {
  Remove-Item $destination
}

Invoke-WebRequest -Uri $source -OutFile $destination
