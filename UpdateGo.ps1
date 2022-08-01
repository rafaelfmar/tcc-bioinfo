$source = 'http://purl.obolibrary.org/obo/go.obo'
$destination = './data/go.obo'

if (Test-Path $destination) {
  Remove-Item $destination
}

Invoke-WebRequest -Uri $source -OutFile $destination
