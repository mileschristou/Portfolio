# Blog Post Framework for mileschristou.co.uk

This document provides the complete structure and formatting guidelines for creating blog posts on the portfolio site.

## Overview

Blog posts are **React JSX components** using **Tailwind CSS** for styling. Each post follows a consistent structure for SEO, readability, and visual design.

## File Structure

- Location: `src/pages/blog/`
- Naming: Kebab-case (e.g., `how-do-soaps-work.jsx`)
- After creating the file, you must:
  1. Import it in `src/App.jsx`
  2. Add a route for it
  3. Add entry to blog listing in `src/pages/Blog.jsx`

## Template Structure

### 1. Imports

```jsx
import { Link } from 'react-router-dom';
import SEO from '../../components/SEO';
```

### 2. Component Setup

```jsx
/*
META DESCRIPTION:
Your meta description here (under 160 characters, includes main keywords)
*/

const YourBlogPost = () => {
  return (
    <>
      {/* SEO and content here */}
    </>
  );
};

export default YourBlogPost;
```

### 3. SEO Component

```jsx
<SEO
  title="Blog Post Title | Miles Christou"
  description="Meta description matching the comment above"
/>
```

### 4. Structured Data - Article Schema

```jsx
<script type="application/ld+json">
  {JSON.stringify({
    "@context": "https://schema.org",
    "@type": "Article",
    "headline": "Full Article Title",
    "description": "Article description",
    "author": {
      "@type": "Person",
      "name": "Miles Christou",
      "jobTitle": "Chemistry Graduate"
    },
    "datePublished": "2026-02-07",
    "dateModified": "2026-02-07",
    "keywords": ["keyword1", "keyword2", "keyword3"],
    "articleSection": "Chemistry",
    "wordCount": 3000
  })}
</script>
```

### 5. Structured Data - FAQ Schema

```jsx
<script type="application/ld+json">
  {JSON.stringify({
    "@context": "https://schema.org",
    "@type": "FAQPage",
    "mainEntity": [
      {
        "@type": "Question",
        "name": "Question here?",
        "acceptedAnswer": {
          "@type": "Answer",
          "text": "Answer here."
        }
      },
      // 4-6 questions total
    ]
  })}
</script>
```

### 6. Article Container

```jsx
<article className="py-12 px-4">
  <div className="max-w-3xl mx-auto">
    {/* Back link */}
    <Link to="/blog" className="inline-block text-gray-600 hover:text-gray-900 mb-8">
      ← Back to Blog
    </Link>

    {/* Header */}
    <header className="mb-12">
      <h1 className="text-4xl md:text-5xl font-bold text-gray-900 mb-4 font-mono leading-tight">
        Your Blog Post Title
      </h1>
      <div className="text-gray-600 text-sm">
        <time dateTime="2026-02-07">February 7, 2026</time>
        <span className="mx-2">•</span>
        <span>12 min read</span>
      </div>
    </header>

    {/* Content */}
    <div className="prose prose-lg max-w-none">
      {/* All content goes here */}
    </div>

    {/* Back to blog link at bottom */}
    <div className="mt-8 pt-8 border-t border-gray-200">
      <Link to="/blog" className="inline-block text-gray-900 font-medium hover:underline">
        ← Back to all posts
      </Link>
    </div>
  </div>
</article>
```

## Content Formatting Guidelines

### Text Elements

**Opening Paragraph (larger text):**
```jsx
<p className="text-xl text-gray-700 leading-relaxed">
  Your engaging opening paragraph here...
</p>
```

**Regular Paragraphs:**
```jsx
<p className="text-gray-700 leading-relaxed mb-4">
  Regular paragraph text...
</p>
```

**Strong Emphasis:**
```jsx
<strong>emphasized text</strong>
```

**Italics:**
```jsx
<em>italicized text</em>
```

### Headings

**H2 (Main Sections):**
```jsx
<h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">
  Section Title
</h2>
```

**H3 (Subsections):**
```jsx
<h3 className="text-2xl font-bold text-gray-900 mt-8 mb-4">
  Subsection Title
</h3>
```

### Special Content Boxes

**Quick Answer Box:**
```jsx
<div className="bg-gray-50 border-l-4 border-gray-900 p-6 my-8">
  <p className="text-gray-700">
    <strong>Quick answer or key concept</strong>: Explanation here...
  </p>
</div>
```

**Key Takeaways Section:**
```jsx
<h3 className="text-2xl font-bold text-gray-900 mt-8 mb-4">Key Takeaways</h3>
<ul className="list-disc pl-6 mb-8 text-gray-700 space-y-2">
  <li>First key point</li>
  <li>Second key point</li>
  <li>Third key point</li>
</ul>
```

**Regular Bullet Lists:**
```jsx
<ul className="list-disc pl-6 mb-4 text-gray-700 space-y-2">
  <li>Item one</li>
  <li>Item two</li>
  <li>Item three</li>
</ul>
```

### Chemical Equations and Formulas

**Simple Equation Box:**
```jsx
<div className="bg-gray-50 border border-gray-200 p-4 my-6 font-mono text-sm">
  <p className="text-center mb-2">[Reactant] + [Reactant] → [Product]</p>
  <p className="text-center text-xs text-gray-600">Equation description</p>
</div>
```

**Larger/More Complex Equation:**
```jsx
<div className="bg-gray-50 border border-gray-200 p-4 my-6 font-mono text-xs text-center">
  <p>Chemical formula or equation here</p>
  <p className="text-gray-600 mt-2">Description or notes</p>
</div>
```

**Chemical Notation:**
- Subscripts: Use Unicode subscripts (₀₁₂₃₄₅₆₇₈₉)
- Superscripts: Use Unicode superscripts (⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻)
- Arrows: Use → for reactions
- Ionic compounds: Use square brackets [NaOH], [CH₃COO⁻ Na⁺]

### Images and Figures

```jsx
<figure className="my-8">
  <img
    src="/image-filename.png"
    alt="Detailed, keyword-rich alt text describing the image content"
    title="Brief title for the image"
    className="w-full max-w-2xl mx-auto"
    loading="lazy"
    width="1000"
    height="800"
  />
  <figcaption className="text-center text-sm text-gray-600 mt-2">
    Caption explaining what the image shows
  </figcaption>
</figure>
```

**Image sizing classes:**
- `max-w-2xl` - Medium images (default)
- `max-w-3xl` - Larger images
- `max-w-4xl` - Very large images

### Links

**Internal Links (to other blog posts):**
```jsx
<Link to="/blog/other-post-slug" className="text-primary-600 hover:text-primary-700 underline">
  link text
</Link>
```

**External Links:**
```jsx
<a href="https://example.com" target="_blank" rel="noopener noreferrer" className="text-primary-600 hover:text-primary-700 underline">
  link text
</a>
```

### FAQ Section (at end of article)

```jsx
<h2 className="text-3xl font-bold text-gray-900 mt-12 mb-4 font-mono">
  Frequently Asked Questions
</h2>

<div className="space-y-6 mb-12">
  <div className="border-l-4 border-gray-400 pl-6">
    <h3 className="text-xl font-bold text-gray-900 mb-2">Question 1?</h3>
    <p className="text-gray-700">
      Answer to question 1...
    </p>
  </div>

  <div className="border-l-4 border-gray-400 pl-6">
    <h3 className="text-xl font-bold text-gray-900 mb-2">Question 2?</h3>
    <p className="text-gray-700">
      Answer to question 2...
    </p>
  </div>

  {/* 4-6 questions total */}
</div>
```

### Author Byline

```jsx
<div className="mt-12 pt-8 border-t border-gray-200">
  <p className="text-gray-600 italic">
    Written by <strong>Miles Christou</strong>, chemistry graduate
  </p>
</div>
```

## Adding the Post to the Site

### 1. Update App.jsx

Add import:
```jsx
import YourBlogPost from './pages/blog/your-blog-post';
```

Add route:
```jsx
<Route path="blog/your-blog-post" element={<YourBlogPost />} />
```

### 2. Update Blog.jsx

Add entry to the `posts` array at the top:
```jsx
{
  title: 'Your Blog Post Title',
  excerpt: 'Brief excerpt (1-2 sentences) describing the post',
  date: 'February 7, 2026',
  readTime: '12 min read',
  slug: 'your-blog-post'
}
```

## Content Guidelines

### Length
- Target: 2,500-3,500 words
- Minimum: 2,000 words for good SEO

### Structure
1. Opening paragraph (engaging, explains what the post covers)
2. Quick Answer section (TL;DR)
3. Key Takeaways (bullet points)
4. Main content sections with H2/H3 headings
5. FAQ section at the end (4-6 questions)
6. Author byline

### Chemistry-Specific
- Use proper chemical notation and Unicode characters
- Bracket ionic compounds in equations
- Link to authoritative scientific sources (journals, CDC, FDA, university sources)
- Include internal links to related chemistry blog posts
- Explain complex concepts clearly for a general audience
- Include diagrams/images where helpful

### SEO Best Practices
- H1 (title) → H2 (sections) → H3 (subsections) hierarchy
- Descriptive alt text on all images (include keywords naturally)
- Meta description under 160 characters
- Target keywords used naturally throughout
- External links to authoritative sources
- Internal links to related content
- Structured data (Article + FAQ schema)
- Date published and modified

### Tone
- Clear and accessible
- Educational but not condescending
- Use active voice
- Explain technical terms when first introduced
- Conversational but professional

## Example File Reference

See `src/pages/blog/how-do-soaps-work.jsx` for a complete working example following all these guidelines.

## Sitemap Updates

After adding a new blog post, update `public/sitemap.xml`:

```xml
<url>
  <loc>https://mileschristou.co.uk/blog/your-blog-post</loc>
  <lastmod>2026-02-07</lastmod>
  <changefreq>monthly</changefreq>
  <priority>0.9</priority>
</url>
```

Then submit the updated sitemap to Google Search Console.

---

**Questions?** Refer to existing blog posts in `src/pages/blog/` for working examples.
