import { DM_Sans, JetBrains_Mono } from 'next/font/google';
import Navbar from '../components/Navbar';
import Footer from '../components/Footer';
import './globals.css';

const dmSans = DM_Sans({
  subsets: ['latin'],
  weight: ['400', '500', '600', '700'],
  style: ['normal', 'italic'],
  variable: '--font-dm-sans',
  display: 'swap',
});

const jetBrainsMono = JetBrains_Mono({
  subsets: ['latin'],
  weight: ['400', '500', '600', '700'],
  variable: '--font-jetbrains-mono',
  display: 'swap',
});

export const metadata = {
  metadataBase: new URL('https://mileschristou.co.uk'),
  title: {
    default: 'Miles Christou | Chemistry Graduate & Analytical Chemist',
    template: '%s | Miles Christou',
  },
  description: 'Portfolio of Miles Christou, a chemistry graduate with expertise in analytical chemistry, GC-MS, HPLC-MS, and PFAS method development.',
  keywords: ['analytical chemistry', 'GC-MS', 'HPLC-MS', 'PFAS', 'scientific writing', 'chemistry graduate'],
  authors: [{ name: 'Miles Christou' }],
  openGraph: {
    type: 'website',
    locale: 'en_GB',
    siteName: 'Miles Christou',
  },
  twitter: {
    card: 'summary_large_image',
  },
};

export default function RootLayout({ children }) {
  return (
    <html lang="en" className={`${dmSans.variable} ${jetBrainsMono.variable}`}>
      <body className="min-h-screen flex flex-col bg-white">
        <Navbar />
        <main className="flex-grow">
          {children}
        </main>
        <Footer />
      </body>
    </html>
  );
}
