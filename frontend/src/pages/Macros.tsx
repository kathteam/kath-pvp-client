import { Fragment, JSX, useEffect } from 'react';
import { AutoAwesome as AutoAwesomeIcon } from '@mui/icons-material';
import { RouteHeader } from '@/components';
import { handleScroll } from '@/utils';

export default function Macros(): JSX.Element {

  // Reset scroll position on initial load
  useEffect(() => {
    handleScroll('Macros');
  }, []);
  
  return (
    <Fragment>
      <RouteHeader
        icon={AutoAwesomeIcon}
        title="Macros"
        description="Record and reuse actions for repetitive tasks to save time and improve accuracy in your workflows."
      />
    </Fragment>
  );
}
