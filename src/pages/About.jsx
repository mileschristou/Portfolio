import SEO from '../components/SEO';

const About = () => {
  const skills = {
    'Laboratory & Analytical Chemistry': [
      'Gas Chromatography-Mass Spectrometry (GC-MS)',
      'High-Performance Liquid Chromatography-Mass Spectrometry (HPLC-MS)',
      'PFAS Method Development',
      'Sample Preparation & Analysis',
      'Quality Control & Data Analysis'
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
        description="Chemistry graduate with experience in analytical chemistry, laboratory work, and scientific communication. Skilled in GC-MS, HPLC-MS, and practical problem-solving."
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
                I am a chemistry graduate with hands-on experience across analytical chemistry,
                method development, and scientific communication. My background spans advanced
                analytical techniques including gas and liquid chromatography-mass spectrometry,
                with a strong foundation in understanding chemical formulations and product chemistry.
              </p>
              <p>
                Throughout my academic journey, I have developed practical skills in laboratory
                work, quality control, and data analysis. During my analytical chemistry modules,
                I worked on PFAS (Per- and Polyfluoroalkyl Substances) method development for
                environmental samples, combining analytical rigor with problem-solving in real-world applications.
              </p>
              <p>
                I combine technical chemistry expertise with clear scientific writing, translating
                complex concepts into accessible content. Beyond traditional lab work, I have an
                interest in web development, creating tools that make scientific concepts more
                approachable.
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
                <h3 className="font-semibold text-primary-900 mb-2">Chemistry & Laboratory Roles</h3>
                <p className="text-primary-800">
                  Analytical chemistry, quality control, product development, formulation,
                  or R&D positions where I can apply practical chemistry skills and
                  problem-solving.
                </p>
              </div>
              <div className="bg-primary-50 border border-primary-200 rounded-lg p-6">
                <h3 className="font-semibold text-primary-900 mb-2">Scientific Communication</h3>
                <p className="text-primary-800">
                  Roles that leverage my ability to translate complex scientific
                  concepts into clear, accessible content for diverse audiences.
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
