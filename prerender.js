import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const distPath = path.join(__dirname, 'dist');
const indexPath = path.join(distPath, 'index.html');

// All routes that need to be pre-rendered
const routes = [
  '/about',
  '/projects',
  '/blog',
  '/contact',
  '/blog/maillard-reaction-cooking',
  '/blog/how-do-soaps-work',
  '/blog/what-are-pfas',
  '/blog/what-are-terpenes',
  '/blog/decaf-coffee-science',
];

// Read the base index.html
const indexHtml = fs.readFileSync(indexPath, 'utf-8');

// Create a directory and index.html for each route
routes.forEach(route => {
  const routePath = path.join(distPath, route);
  const routeIndexPath = path.join(routePath, 'index.html');

  // Create directory if it doesn't exist
  fs.mkdirSync(routePath, { recursive: true });

  // Copy index.html to the route directory
  fs.writeFileSync(routeIndexPath, indexHtml);

  console.log(`✓ Pre-rendered: ${route}`);
});

console.log('\n✓ Pre-rendering complete!');
