import Link from 'next/link';

export default function NotFound() {
  return (
    <div className="py-20 px-4">
      <div className="max-w-4xl mx-auto text-center">
        <h1 className="text-6xl font-bold text-gray-900 mb-4 font-mono">404</h1>
        <p className="text-xl text-gray-700 mb-8">
          This page doesn&apos;t exist.
        </p>
        <Link
          href="/"
          className="inline-block px-8 py-3 bg-gray-900 text-white font-medium hover:bg-gray-800 transition-colors"
        >
          Go Home
        </Link>
      </div>
    </div>
  );
}
