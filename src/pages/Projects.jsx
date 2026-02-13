import SEO from '../components/SEO';

const Projects = () => {
  const projects = [
    {
      title: 'SwapTheStyle (HomeMockups)',
      description: 'AI-powered interior and exterior design tool that transforms any space instantly with 58+ design styles. Upload photos of kitchens, living rooms, bedrooms, bathrooms, gardens, or patios and receive professional AI-generated redesigns in seconds. Built with Next.js and React, integrating advanced AI processing for instant design generation across multiple aesthetic categories.',
      technologies: ['Next.js', 'React', 'AI', 'Image Processing', 'Design'],
      link: 'https://homemockups.com'
    },
    {
      title: 'Ironman Meta',
      description: 'Practical, opinionated guides for Old School RuneScape Ironman mode players. Features stage-by-stage progression guides from Tutorial Island through end-game, ironman-specific skilling methods, supply acquisition strategies, and in-depth coverage of specific regions and activities. Provides clear, actionable advice grounded in practical experience rather than wiki-style data dumps.',
      technologies: ['HTML', 'CSS', 'JavaScript', 'Content Management'],
      github: 'https://github.com/mileschristou/ironman-meta',
      link: 'https://ironmanmeta.com'
    },
    {
      title: 'Rust Server Tracker',
      description: 'Full-stack web application for tracking Rust game server wipes and statistics. Features real-time server data from BattleMetrics API, AI-powered server search using Anthropic, advanced filtering by player count, region, and server type, and automated wipe detection with historical tracking. Built with Node.js/Express backend, React frontend, and PostgreSQL database with scheduled data scraping.',
      technologies: ['React', 'Node.js', 'Express', 'PostgreSQL', 'AI', 'BattleMetrics API'],
      github: 'https://github.com/mileschristou/rust-server-tracker',
      link: '#'
    },
    {
      title: 'NMR Peak Directory',
      description: 'A comprehensive searchable database of AI-generated annotated NMR spectra. Features include compound search, filterable by nucleus type (1H, 13C), solvent, and functional groups. Each entry displays labeled chemical shift peaks with detailed descriptions of characteristic patterns, designed to help chemistry students and researchers quickly identify NMR spectral features.',
      technologies: ['React', 'JavaScript', 'AI', 'NMR Spectroscopy', 'Data Visualization'],
      github: 'https://github.com/mileschristou/NMR-Peaks-Directory',
      link: 'https://nmr-peaks-directory-5ojd0xvex-miles-projects-89a76c0f.vercel.app/'
    },
    {
      title: 'AI d-Electron Counter',
      description: 'An intelligent web application that calculates d-electron counts for transition metal complexes using AI. Users can input a metal with its oxidation state or a complete coordination complex formula, and the AI provides the d-electron count with detailed step-by-step reasoning. Features both a calculator mode and learning resources to help students master this fundamental concept in inorganic chemistry.',
      technologies: ['React', 'AI', 'Inorganic Chemistry', 'Authentication', 'Educational'],
      link: 'https://d-electron-counter.vercel.app/'
    }
  ];

  return (
    <>
      <SEO
        title="Projects | Miles Christou"
        description="Explore Miles Christou's web applications including AI-powered design tools, gaming guides, and chemistry educational resources - from SwapTheStyle interior design to Ironman Meta and NMR spectroscopy tools."
      />

      <div className="py-20 px-4">
        <div className="max-w-4xl mx-auto">
          {/* Header */}
          <div className="mb-16">
            <h1 className="text-4xl md:text-5xl font-bold text-gray-900 mb-4 font-mono">
              Projects
            </h1>
            <p className="text-lg text-gray-700 max-w-2xl">
              A collection of web applications spanning AI-powered design tools, gaming resources, and chemistry education.
            </p>
          </div>

          {/* Projects */}
          <div className="space-y-16">
            {projects.map((project, index) => (
              <div key={index} className="border-l-4 border-gray-900 pl-8 py-2">
                <h2 className="text-2xl font-bold text-gray-900 mb-3 font-mono">
                  {project.title}
                </h2>
                <p className="text-gray-700 mb-4 leading-relaxed">
                  {project.description}
                </p>
                <div className="flex flex-wrap gap-2 mb-4">
                  {project.technologies.map((tech, idx) => (
                    <span
                      key={idx}
                      className="text-sm text-gray-600 bg-gray-100 px-3 py-1 font-mono"
                    >
                      {tech}
                    </span>
                  ))}
                </div>
                <div className="flex gap-4">
                  {project.link && project.link !== '#' && (
                    <a
                      href={project.link}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="inline-block text-gray-900 font-medium hover:underline"
                    >
                      View Live Demo →
                    </a>
                  )}
                  {project.github && (
                    <a
                      href={project.github}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="inline-block text-gray-900 font-medium hover:underline"
                    >
                      View on GitHub →
                    </a>
                  )}
                  {(!project.link || project.link === '#') && !project.github && (
                    <span className="inline-block text-gray-500 font-medium">
                      Coming Soon
                    </span>
                  )}
                </div>
              </div>
            ))}
          </div>

          {/* Note */}
          <div className="mt-16 p-6 bg-gray-50 border border-gray-200">
            <p className="text-gray-700">
              More projects in development. Check back soon for updates.
            </p>
          </div>
        </div>
      </div>
    </>
  );
};

export default Projects;
