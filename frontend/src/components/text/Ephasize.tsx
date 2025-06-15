import { JSX } from 'react'
import { Typography } from '@mui/material'
import { useTheme } from '@mui/material/styles'

export default function Ephasize(props: { label: string, text: string }) : JSX.Element {
  const theme = useTheme()

  return (
    <>
      <Typography variant="body2" sx={{ fontWeight: 600, color: theme.palette.text.primary, mr: 1 }}>
        {props.label}
      </Typography>
      <Typography variant="body2" sx={{ color: theme.palette.text.primary, mr: 1 }}>
        {props.text}
      </Typography>
    </>
  )
}
