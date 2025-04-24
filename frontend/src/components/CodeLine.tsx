import { Box, Typography, useTheme } from '@mui/material';
import { useEffect, useRef, useState } from 'react';

interface CodeLineProps {
  lineNumber: number;
  content: string;
  bgcolor?: string;
}

export default function CodeLine(props: CodeLineProps) {
  const { lineNumber, content, bgcolor } = props;
  const theme = useTheme();
  const contentRef = useRef<HTMLDivElement>(null);
  const [isWrapped, setIsWrapped] = useState(false);

  useEffect(() => {
    const el = contentRef.current;
    if (el) {
      const lineHeight = parseFloat(getComputedStyle(el).lineHeight || '20');
      const actualHeight = el.getBoundingClientRect().height;
      setIsWrapped(actualHeight > lineHeight + 1); // 1px buffer
    }
  }, [content]);

  return (
    <Box
      sx={{
        display: 'flex',
        alignItems: 'flex-start',
        width: '100%',
        bgcolor: bgcolor || theme.palette.background.default,
      }}
    >
      {/* Line Number */}
      <Typography
        sx={{
          width: 40,
          minWidth: 40,
          color: theme.palette.text.secondary,
          textAlign: 'right',
          pr: 1,
          userSelect: 'none',
          flexShrink: 0,
        }}
      >
        {lineNumber}
      </Typography>

      {/* Line Content */}
      <Box
        ref={contentRef}
        sx={{
          flex: 1,
          whiteSpace: 'pre-wrap',
          wordBreak: 'break-word',
          color: theme.palette.text.primary,
          position: 'relative',
          fontFamily: theme.typography.fontFamily,
          fontSize: theme.typography.body1.fontSize
        }}
      >
        {content}
        {isWrapped && (
          <Box
            component="span"
            sx={{
              position: 'absolute',
              bottom: 0,
              right: 0,
              px: 0.5,
              fontWeight: 'bold',
              backgroundColor: 'transparent',
              color: theme.palette.mode === 'light' ? theme.palette.info.main : theme.palette.info.dark,
              fontSize: '0.75em',
            }}
          >
            ↩
          </Box>
        )}
      </Box>
    </Box>
  );
}
