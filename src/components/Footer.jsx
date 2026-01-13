const Footer = () => {
  const currentYear = new Date().getFullYear();

  return (
    <footer className="bg-white border-t border-gray-200 mt-auto">
      <div className="max-w-6xl mx-auto px-4 sm:px-6 lg:px-8 py-12">
        <div className="flex flex-col md:flex-row justify-between items-center space-y-4 md:space-y-0">
          <div className="text-center md:text-left">
            <p className="text-sm text-gray-600">
              &copy; {currentYear} Miles Christou
            </p>
          </div>

          <div className="flex gap-6">
            <a href="/about" className="text-sm text-gray-600 hover:text-gray-900 transition-colors">
              About
            </a>
            <a href="/projects" className="text-sm text-gray-600 hover:text-gray-900 transition-colors">
              Projects
            </a>
            <a
              href="https://www.linkedin.com/in/miles-christou-742909333/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-sm text-gray-600 hover:text-gray-900 transition-colors"
            >
              LinkedIn
            </a>
          </div>
        </div>
      </div>
    </footer>
  );
};

export default Footer;
