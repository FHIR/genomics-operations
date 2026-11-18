import { useEffect, useState } from 'react';

/**
 * Height of the sticky patient banner, kept current as the window resizes.
 */
export function useBannerHeight() {
  const [height, setHeight] = useState(0);
  useEffect(() => {
    const banner = document.getElementById('patient-banner');
    if (!banner) {
      return;
    }
    const update = () => setHeight(banner.offsetHeight);
    update();
    const observer = new ResizeObserver(update);
    observer.observe(banner);
    return () => observer.disconnect();
  }, []);
  return height;
}
