# Miles Christou - Portfolio Website

🌐 **Live Site:** https://portfolio-smoky-kappa-53.vercel.app/

A professional portfolio website built with React, Vite, Tailwind CSS, and React Router v6.

## Features

- Clean, professional design with chemistry-inspired color palette
- Fully responsive (mobile, tablet, desktop)
- SEO-optimized with custom meta tags and structured data (FAQ schema, image schema)
- Technical blog posts with scientific citations and molecular diagrams
- Project showcase featuring deployed web applications
- Fast build and development with Vite
- Deployed on Vercel

## Project Structure

```
portfolio/
├── src/
│   ├── components/
│   │   ├── Navbar.jsx       # Navigation with mobile menu
│   │   ├── Footer.jsx       # Footer with social links
│   │   ├── Layout.jsx       # Main layout wrapper
│   │   └── SEO.jsx          # Custom SEO component
│   ├── pages/
│   │   ├── Home.jsx         # Landing page
│   │   ├── About.jsx        # Background & skills
│   │   ├── Projects.jsx     # Project showcase
│   │   ├── Blog.jsx         # Blog listing page
│   │   ├── Contact.jsx      # Contact information
│   │   └── blog/            # Individual blog posts
│   │       ├── decaf-coffee-science.jsx
│   │       └── what-are-terpenes.jsx
│   ├── App.jsx              # React Router setup
│   ├── main.jsx
│   └── index.css            # Tailwind configuration
├── vercel.json              # Vercel SPA routing config
└── package.json
```

## Getting Started

### Development

Start the development server:
```bash
cd ~/Desktop/portfolio
npm run dev
```

The site will be available at `http://localhost:5173/`

### Build for Production

```bash
npm run build
```

The built files will be in the `dist/` folder.

### Preview Production Build

```bash
npm run preview
```

## Deployment to Vercel

1. Push your code to a GitHub repository
2. Connect your repository to Vercel
3. Vercel will automatically detect Vite settings
4. Deploy!

Or use the Vercel CLI:
```bash
npm install -g vercel
vercel
```

## Customization

### Update Personal Information

1. **LinkedIn URL**: Update in `src/components/Footer.jsx` and `src/pages/Contact.jsx`
2. **Email**: Update in `src/pages/Contact.jsx`
3. **Project Links**: Update placeholder links in `src/pages/Projects.jsx`

### Add Blog Posts

The blog structure is ready. To add blog posts:
1. Create blog post components in `src/pages/blog/`
2. Update `src/pages/Blog.jsx` to list and link to posts
3. Each post can use the `<SEO>` component for optimization

### Customize Colors

Edit the color palette in `src/index.css` under the `@theme` section.

## SEO

The custom SEO component (`src/components/SEO.jsx`) handles:
- Dynamic page titles
- Meta descriptions
- Open Graph tags for social sharing
- Twitter Card tags

Update SEO props in each page component to customize meta tags.

## Technologies

- **React 19** - UI framework
- **Vite 7** - Build tool and dev server
- **Tailwind CSS 4** - Utility-first CSS framework
- **React Router v6** - Client-side routing
- **Custom SEO component** - React 19 compatible meta tag management

## License

Personal portfolio - all rights reserved.
