import { BrowserRouter, Routes, Route } from 'react-router-dom';
import Layout from './components/Layout';
import Home from './pages/Home';
import About from './pages/About';
import Projects from './pages/Projects';
import Blog from './pages/Blog';
import Contact from './pages/Contact';
import DecafCoffeeScience from './pages/blog/decaf-coffee-science';
import WhatAreTerpenes from './pages/blog/what-are-terpenes';
import WhatArePFAS from './pages/blog/what-are-pfas';
import HowDoSoapsWork from './pages/blog/how-do-soaps-work';

function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route path="/" element={<Layout />}>
          <Route index element={<Home />} />
          <Route path="about" element={<About />} />
          <Route path="projects" element={<Projects />} />
          <Route path="blog" element={<Blog />} />
          <Route path="blog/decaf-coffee-science" element={<DecafCoffeeScience />} />
          <Route path="blog/what-are-terpenes" element={<WhatAreTerpenes />} />
          <Route path="blog/what-are-pfas" element={<WhatArePFAS />} />
          <Route path="blog/how-do-soaps-work" element={<HowDoSoapsWork />} />
          <Route path="contact" element={<Contact />} />
        </Route>
      </Routes>
    </BrowserRouter>
  );
}

export default App;
