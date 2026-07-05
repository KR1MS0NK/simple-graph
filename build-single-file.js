const fs = require('fs');
const path = require('path');

console.log("Starting build process for standalone HTML...");

try {
  const indexHtmlPath = path.join(__dirname, 'index.html');
  const stylesCssPath = path.join(__dirname, 'styles.css');
  const appJsPath = path.join(__dirname, 'app.js');
  const outputPath = path.join(__dirname, 'dist', 'Graph-Analysis-Standalone.html');

  // Read files
  let indexHtml = fs.readFileSync(indexHtmlPath, 'utf8');
  const stylesCss = fs.readFileSync(stylesCssPath, 'utf8');
  const appJs = fs.readFileSync(appJsPath, 'utf8');

  // 1. Embed CSS
  const cssPlaceholder = '<link rel="stylesheet" href="styles.css">';
  if (indexHtml.includes(cssPlaceholder)) {
    indexHtml = indexHtml.replace(cssPlaceholder, `<style>\n${stylesCss}\n</style>`);
    console.log("- Embedded styles.css successfully.");
  } else {
    console.warn("! Warning: styles.css link placeholder not found in index.html.");
  }

  // 2. Embed JS
  const jsPlaceholder = '<script src="app.js"></script>';
  if (indexHtml.includes(jsPlaceholder)) {
    indexHtml = indexHtml.replace(jsPlaceholder, `<script>\n${appJs}\n</script>`);
    console.log("- Embedded app.js successfully.");
  } else {
    console.warn("! Warning: app.js script placeholder not found in index.html.");
  }

  // Ensure output directory exists
  const distDir = path.join(__dirname, 'dist');
  if (!fs.existsSync(distDir)) {
    fs.mkdirSync(distDir);
  }

  // Write standalone single-file
  fs.writeFileSync(outputPath, indexHtml, 'utf8');
  console.log(`\nSUCCESS! Standalone HTML application generated at:\n${outputPath}`);
  process.exit(0);

} catch (err) {
  console.error("ERROR generating standalone HTML:", err);
  process.exit(1);
}
