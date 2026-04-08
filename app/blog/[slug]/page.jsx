import Link from 'next/link';
import { notFound } from 'next/navigation';
import { posts, getPostBySlug, getAllSlugs } from '../../../lib/posts';

import DecafCoffeeScienceContent from '../../../components/blog/decaf-coffee-science';
import WhatAreTerpenesContent from '../../../components/blog/what-are-terpenes';
import WhatArePFASContent from '../../../components/blog/what-are-pfas';
import HowDoSoapsWorkContent from '../../../components/blog/how-do-soaps-work';
import MaillardReactionCookingContent from '../../../components/blog/maillard-reaction-cooking';
import AnalyticalMethodsContent from '../../../components/blog/analytical-methods-icp-ic-gc';

const contentComponents = {
  'decaf-coffee-science': DecafCoffeeScienceContent,
  'what-are-terpenes': WhatAreTerpenesContent,
  'what-are-pfas': WhatArePFASContent,
  'how-do-soaps-work': HowDoSoapsWorkContent,
  'maillard-reaction-cooking': MaillardReactionCookingContent,
  'analytical-methods-icp-ic-gc': AnalyticalMethodsContent,
};

export function generateStaticParams() {
  return getAllSlugs().map((slug) => ({ slug }));
}

export async function generateMetadata({ params }) {
  const { slug } = await params;
  const post = getPostBySlug(slug);
  if (!post) return {};

  return {
    title: post.title,
    description: post.description,
    keywords: post.keywords,
    openGraph: {
      title: post.title,
      description: post.description,
      type: 'article',
      publishedTime: post.dateISO,
      modifiedTime: post.dateModified,
      authors: ['Miles Christou'],
    },
    alternates: {
      canonical: `/blog/${post.slug}`,
    },
  };
}

export default async function BlogPost({ params }) {
  const { slug } = await params;
  const post = getPostBySlug(slug);

  if (!post) {
    notFound();
  }

  const ContentComponent = contentComponents[slug];

  if (!ContentComponent) {
    notFound();
  }

  const articleSchema = {
    '@context': 'https://schema.org',
    '@type': 'Article',
    headline: post.title,
    description: post.description,
    author: {
      '@type': 'Person',
      name: 'Miles Christou',
      jobTitle: 'Chemistry Graduate',
    },
    datePublished: post.dateISO,
    dateModified: post.dateModified,
    keywords: post.keywords,
    articleSection: post.articleSection,
    wordCount: post.wordCount,
    url: `https://mileschristou.co.uk/blog/${post.slug}`,
    publisher: {
      '@type': 'Person',
      name: 'Miles Christou',
    },
    ...(post.imageSchemas ? { image: post.imageSchemas } : {}),
  };

  const faqSchema = post.faqSchema ? {
    '@context': 'https://schema.org',
    '@type': 'FAQPage',
    mainEntity: post.faqSchema.map((faq) => ({
      '@type': 'Question',
      name: faq.question,
      acceptedAnswer: {
        '@type': 'Answer',
        text: faq.answer,
      },
    })),
  } : null;

  return (
    <>
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(articleSchema) }}
      />
      {faqSchema && (
        <script
          type="application/ld+json"
          dangerouslySetInnerHTML={{ __html: JSON.stringify(faqSchema) }}
        />
      )}

      <article className="py-12 px-4">
        <div className="max-w-3xl mx-auto">
          {/* Back link */}
          <Link href="/blog" className="inline-block text-gray-600 hover:text-gray-900 mb-8">
            ← Back to Blog
          </Link>

          {/* Article header */}
          <header className="mb-12">
            <h1 className="text-4xl md:text-5xl font-bold text-gray-900 mb-4 font-mono leading-tight">
              {post.title}
            </h1>
            <div className="text-gray-600 text-sm">
              <time dateTime={post.dateISO}>{post.date}</time>
              <span className="mx-2">•</span>
              <span>{post.readTime}</span>
            </div>
          </header>

          {/* Article content */}
          <ContentComponent />

          {/* Back to blog link */}
          <div className="mt-8 pt-8 border-t border-gray-200">
            <Link href="/blog" className="inline-block text-gray-900 font-medium hover:underline">
              ← Back to all posts
            </Link>
          </div>
        </div>
      </article>
    </>
  );
}
