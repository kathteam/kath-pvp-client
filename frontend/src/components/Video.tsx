import { JSX } from 'react';
import { Column, MediaContainer, Row } from '@/components/core';

export interface VideoProps {
  src: string;
}

export default function Video({ src }: VideoProps): JSX.Element {
  return (
    <Column>
      <MediaContainer>
        <Row sx={{ p: 0 }}>
          <video
            src={src}
            controls
            style={{ width: '100%', height: 'auto', borderRadius: '12px' }}
          />   
        </Row>
      </MediaContainer>
    </Column>
  );
}
