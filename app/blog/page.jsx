import Link from 'next/link';
import { posts } from '../../lib/posts';

export const metadata = {
  title: 'Blog',
  description: 'Scientific writing samples and articles by Miles Christou covering analytical chemistry, method development, and scientific communication.',
};

export default function Blog() {
  return (
    <div className="py-20 px-4">
      <div className="max-w-4xl mx-auto">
        {/* Header */}
        <div className="mb-16">
          <h1 className="text-4xl md:text-5xl font-bold text-gray-900 mb-4 font-mono">
            Blog
          </h1>
          <p className="text-lg text-gray-700">
            Scientific writing samples, insights on analytical chemistry, and explorations of method development and research communication.
          </p>
        </div>

        {/* Blog posts */}
        <div className="space-y-12">
          {posts.map((post, index) => (
            <article key={index} className="border-l-4 border-gray-900 pl-8 py-2">
              <Link href={`/blog/${post.slug}`} className="group">
                <h2 className="text-2xl md:text-3xl font-bold text-gray-900 mb-3 font-mono group-hover:text-gray-600 transition-colors">
                  {post.title}
                </h2>
                <div className="text-sm text-gray-600 mb-3">
                  {post.date} <span className="mx-2">•</span> {post.readTime}
                </div>
                <p className="text-gray-700 leading-relaxed mb-4">
                  {post.excerpt}
                </p>
                <span className="inline-block text-gray-900 font-medium group-hover:underline">
                  Read more →
                </span>
              </Link>
            </article>
          ))}
        </div>
      </div>
    </div>
  );
}
