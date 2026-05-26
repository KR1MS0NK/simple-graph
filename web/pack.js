const packager = require('electron-packager');

console.log("Starting programmatic packager...");
packager({
  dir: '.',
  name: 'Graph-Analysis',
  platform: 'win32',
  arch: 'x64',
  out: 'dist',
  overwrite: true
}).then(appPaths => {
  console.log("SUCCESS! Wrote app to:", appPaths);
  process.exit(0);
}).catch(err => {
  console.error("ERROR during packaging:", err);
  process.exit(1);
});
