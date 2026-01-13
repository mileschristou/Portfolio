import SEO from '../components/SEO';

const About = () => {
  const skills = {
    'Analytical Chemistry': [
      'Gas Chromatography-Mass Spectrometry (GC-MS)',
      'High-Performance Liquid Chromatography-Mass Spectrometry (HPLC-MS)',
      'PFAS Method Development',
      'Sample Preparation & Analysis',
      'Quality Control & Validation'
    ],
    'Scientific Writing': [
      'Technical Documentation',
      'Research Communication',
      'Method Development Reports',
      'Scientific Content Creation'
    ],
    'Web Development': [
      'React & JavaScript',
      'Educational Tool Development',
      'Data Visualization',
      'Responsive Web Design'
    ]
  };

  return (
    <>
      <SEO
        title="About | Miles Christou"
        description="Learn about Miles Christou's background in analytical chemistry, expertise in GC-MS and HPLC-MS, and passion for scientific writing."
      />

      <div className="py-16 px-4">
        <div className="max-w-4xl mx-auto">
          {/* Header */}
          <div className="mb-12">
            <h1 className="text-4xl font-bold text-gray-900 mb-4">About Me</h1>
            <div className="w-20 h-1 bg-primary-600 mb-8"></div>
          </div>

          {/* Background */}
          <section className="mb-12">
            <h2 className="text-2xl font-semibold text-gray-900 mb-4">Background</h2>
            <div className="prose prose-lg text-gray-700 space-y-4">
              <p>
                I am a chemistry graduate with a strong foundation in analytical chemistry
                and a passion for scientific communication. My expertise lies in advanced
                analytical techniques, particularly in gas and liquid chromatography coupled
                with mass spectrometry.
              </p>
              <p>
                Throughout my academic and professional journey, I have developed specialized
                skills in PFAS (Per- and Polyfluoroalkyl Substances) method development,
                contributing to important environmental analysis work. I combine technical
                expertise with clear scientific writing to communicate complex findings
                effectively.
              </p>
              <p>
                Beyond traditional chemistry work, I have cultivated an interest in web
                development, creating educational tools that bridge the gap between complex
                scientific concepts and accessible learning resources.
              </p>
            </div>
          </section>

          {/* Skills */}
          <section className="mb-12">
            <h2 className="text-2xl font-semibold text-gray-900 mb-6">Skills & Expertise</h2>
            <div className="space-y-6">
              {Object.entries(skills).map(([category, skillList]) => (
                <div key={category} className="bg-gray-50 rounded-lg p-6">
                  <h3 className="text-xl font-semibold text-gray-900 mb-4">{category}</h3>
                  <div className="flex flex-wrap gap-2">
                    {skillList.map((skill, index) => (
                      <span
                        key={index}
                        className="bg-white text-gray-700 px-4 py-2 rounded-md border border-gray-200 text-sm"
                      >
                        {skill}
                      </span>
                    ))}
                  </div>
                </div>
              ))}
            </div>
          </section>

          {/* Career Goals */}
          <section>
            <h2 className="text-2xl font-semibold text-gray-900 mb-4">Career Interests</h2>
            <p className="text-gray-700 text-lg mb-4">
              I am actively seeking opportunities in:
            </p>
            <div className="grid md:grid-cols-2 gap-4">
              <div className="bg-primary-50 border border-primary-200 rounded-lg p-6">
                <h3 className="font-semibold text-primary-900 mb-2">Analytical Chemistry Roles</h3>
                <p className="text-primary-800">
                  Laboratory positions focusing on advanced analytical techniques,
                  method development, and environmental analysis.
                </p>
              </div>
              <div className="bg-primary-50 border border-primary-200 rounded-lg p-6">
                <h3 className="font-semibold text-primary-900 mb-2">Scientific Writing Positions</h3>
                <p className="text-primary-800">
                  Roles that leverage my ability to communicate complex scientific
                  concepts clearly and effectively to diverse audiences.
                </p>
              </div>
            </div>
          </section>
        </div>
      </div>
    </>
  );
};

export default About;
