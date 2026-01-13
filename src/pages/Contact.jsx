import SEO from '../components/SEO';

const Contact = () => {
  return (
    <>
      <SEO
        title="Contact | Miles Christou"
        description="Get in touch with Miles Christou for analytical chemistry opportunities, scientific writing projects, or collaboration inquiries."
      />

      <div className="py-16 px-4">
        <div className="max-w-4xl mx-auto">
          {/* Header */}
          <div className="mb-12">
            <h1 className="text-4xl font-bold text-gray-900 mb-4">Contact</h1>
            <div className="w-20 h-1 bg-primary-600 mb-6"></div>
            <p className="text-lg text-gray-700">
              I'm interested in analytical chemistry roles, scientific writing positions,
              and collaborative projects. Let's connect!
            </p>
          </div>

          {/* Contact Cards */}
          <div className="grid md:grid-cols-2 gap-8 mb-12">
            {/* LinkedIn */}
            <div className="bg-gradient-to-br from-primary-50 to-white border border-primary-200 rounded-lg p-8">
              <div className="flex items-center mb-4">
                <div className="bg-primary-600 p-3 rounded-lg">
                  <svg className="w-6 h-6 text-white" fill="currentColor" viewBox="0 0 24 24">
                    <path d="M19 0h-14c-2.761 0-5 2.239-5 5v14c0 2.761 2.239 5 5 5h14c2.762 0 5-2.239 5-5v-14c0-2.761-2.238-5-5-5zm-11 19h-3v-11h3v11zm-1.5-12.268c-.966 0-1.75-.79-1.75-1.764s.784-1.764 1.75-1.764 1.75.79 1.75 1.764-.783 1.764-1.75 1.764zm13.5 12.268h-3v-5.604c0-3.368-4-3.113-4 0v5.604h-3v-11h3v1.765c1.396-2.586 7-2.777 7 2.476v6.759z" />
                  </svg>
                </div>
                <h2 className="text-2xl font-semibold text-gray-900 ml-4">LinkedIn</h2>
              </div>
              <p className="text-gray-700 mb-4">
                Connect with me on LinkedIn to view my professional profile and experience.
              </p>
              <a
                href="https://www.linkedin.com/in/miles-christou-742909333/"
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center bg-primary-600 text-white px-6 py-3 rounded-lg font-medium hover:bg-primary-700 transition-colors"
              >
                View LinkedIn Profile
                <svg className="w-4 h-4 ml-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
                </svg>
              </a>
            </div>

            {/* Email */}
            <div className="bg-gradient-to-br from-gray-50 to-white border border-gray-200 rounded-lg p-8">
              <div className="flex items-center mb-4">
                <div className="bg-gray-700 p-3 rounded-lg">
                  <svg className="w-6 h-6 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 8l7.89 5.26a2 2 0 002.22 0L21 8M5 19h14a2 2 0 002-2V7a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" />
                  </svg>
                </div>
                <h2 className="text-2xl font-semibold text-gray-900 ml-4">Email</h2>
              </div>
              <p className="text-gray-700 mb-4">
                Reach out directly for opportunities or inquiries.
              </p>
              <a
                href="mailto:mileschristou@gmail.com"
                className="inline-flex items-center bg-gray-700 text-white px-6 py-3 rounded-lg font-medium hover:bg-gray-800 transition-colors"
              >
                Send an Email
                <svg className="w-4 h-4 ml-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M14 5l7 7m0 0l-7 7m7-7H3" />
                </svg>
              </a>
            </div>
          </div>

          {/* Interests Section */}
          <div className="bg-white border border-gray-200 rounded-lg p-8">
            <h2 className="text-2xl font-semibold text-gray-900 mb-4">
              Open to Opportunities In:
            </h2>
            <div className="space-y-3">
              {[
                'Analytical Chemistry Positions (GC-MS, HPLC-MS, PFAS Analysis)',
                'Scientific Writing & Technical Communication',
                'Method Development & Validation Projects',
                'Research & Development Roles',
                'Educational Content Creation'
              ].map((opportunity, index) => (
                <div key={index} className="flex items-start">
                  <svg
                    className="w-6 h-6 text-primary-600 mr-3 flex-shrink-0 mt-0.5"
                    fill="none"
                    stroke="currentColor"
                    viewBox="0 0 24 24"
                  >
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth={2}
                      d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z"
                    />
                  </svg>
                  <span className="text-gray-700 text-lg">{opportunity}</span>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </>
  );
};

export default Contact;
