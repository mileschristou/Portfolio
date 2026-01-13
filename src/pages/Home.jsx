import { Link } from 'react-router-dom';
import SEO from '../components/SEO';

const Home = () => {
  return (
    <>
      <SEO
        title="Miles Christou | Chemistry Graduate & Analytical Chemist"
        description="Portfolio of Miles Christou, a First Class Chemistry graduate from Newcastle University combining analytical chemistry with web development."
      />

      <div className="min-h-[calc(100vh-4rem)]">
        {/* Hero Section - minimal */}
        <section className="py-20 md:py-28 px-4 bg-gray-50">
          <div className="max-w-4xl mx-auto">
            <div className="mb-6">
              <span className="text-sm text-gray-600 font-mono">
                Newcastle University
              </span>
            </div>

            <h1 className="text-5xl md:text-7xl font-bold text-gray-900 mb-6 font-mono">
              Miles Christou
            </h1>

            <div className="space-y-3 mb-8">
              <p className="text-xl md:text-2xl text-gray-700">
                First Class BSc Chemistry Graduate
              </p>
              <p className="text-lg md:text-xl text-gray-600">
                Building Digital Projects
              </p>
            </div>

            <p className="text-lg text-gray-700 mb-10 max-w-2xl leading-relaxed">
              I'm Miles, a First Class Chemistry graduate from Newcastle University with
              a passion for practical and analytical chemistry. I enjoy building digital
              projects that allow me to combine my scientific background with web development.
            </p>

            <div className="flex flex-col sm:flex-row gap-4">
              <Link
                to="/projects"
                className="inline-block px-8 py-3 bg-gray-900 text-white font-medium hover:bg-gray-800 transition-colors"
              >
                View Projects →
              </Link>
              <Link
                to="/about"
                className="inline-block px-8 py-3 border-2 border-gray-900 text-gray-900 font-medium hover:bg-gray-50 transition-colors"
              >
                About Me
              </Link>
            </div>
          </div>
        </section>

        {/* What I Do Section */}
        <section className="py-20 px-4 bg-white">
          <div className="max-w-4xl mx-auto">
            <h2 className="text-3xl md:text-4xl font-bold text-gray-900 mb-12 font-mono">
              What I Do
            </h2>

            <div className="space-y-12">
              {/* Analytical Chemistry */}
              <div className="border-l-4 border-gray-900 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-3 font-mono">
                  Analytical Chemistry
                </h3>
                <p className="text-gray-700 leading-relaxed">
                  Training in advanced analytical techniques and method development for precision chemical analysis.
                </p>
              </div>

              {/* Technical Communication */}
              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-3 font-mono">
                  Technical Communication
                </h3>
                <p className="text-gray-700 leading-relaxed">
                  Clear communication of complex scientific concepts through writing and digital content.
                </p>
              </div>

              {/* Web Development */}
              <div className="border-l-4 border-gray-400 pl-6">
                <h3 className="text-xl font-bold text-gray-900 mb-3 font-mono">
                  Web Development
                </h3>
                <p className="text-gray-700 leading-relaxed">
                  Building educational tools and applications that bridge chemistry and technology.
                </p>
              </div>
            </div>
          </div>
        </section>

        {/* Simple CTA */}
        <section className="py-16 px-4 bg-gray-900 text-white">
          <div className="max-w-4xl mx-auto text-center">
            <h2 className="text-2xl md:text-3xl font-bold mb-4 font-mono">
              Interested in working together?
            </h2>
            <p className="text-lg text-gray-300 mb-8">
              Open to analytical chemistry roles and collaborative projects.
            </p>
            <Link
              to="/contact"
              className="inline-block px-8 py-3 bg-white text-gray-900 font-medium hover:bg-gray-100 transition-colors"
            >
              Get in Touch
            </Link>
          </div>
        </section>
      </div>
    </>
  );
};

export default Home;
