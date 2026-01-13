import { useEffect } from 'react';

const SEO = ({
  title = 'Miles Christou | Chemistry Graduate & Analytical Chemist',
  description = 'Portfolio of Miles Christou, a chemistry graduate with expertise in analytical chemistry, GC-MS, HPLC-MS, and PFAS method development.',
  keywords = 'analytical chemistry, GC-MS, HPLC-MS, PFAS, scientific writing, chemistry graduate',
  type = 'website',
  url = window.location.href
}) => {
  useEffect(() => {
    // Update document title
    document.title = title;

    // Update or create meta tags
    const updateMetaTag = (name, content, attribute = 'name') => {
      let element = document.querySelector(`meta[${attribute}="${name}"]`);
      if (!element) {
        element = document.createElement('meta');
        element.setAttribute(attribute, name);
        document.head.appendChild(element);
      }
      element.setAttribute('content', content);
    };

    // Standard meta tags
    updateMetaTag('description', description);
    updateMetaTag('keywords', keywords);

    // Open Graph tags
    updateMetaTag('og:title', title, 'property');
    updateMetaTag('og:description', description, 'property');
    updateMetaTag('og:type', type, 'property');
    updateMetaTag('og:url', url, 'property');

    // Twitter Card tags
    updateMetaTag('twitter:card', 'summary_large_image');
    updateMetaTag('twitter:title', title);
    updateMetaTag('twitter:description', description);
  }, [title, description, keywords, type, url]);

  return null;
};

export default SEO;
