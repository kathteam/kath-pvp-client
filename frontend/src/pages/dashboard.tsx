import { JSX } from "react"
import { useNavigate } from "react-router-dom"

export default function Dashboard(): JSX.Element {
  const navigate = useNavigate()

  return (
    <div>
      <h1>Dashboard Page</h1>
      <p>This is the dashboard page of our application.</p>
      <button onClick={() => window.pywebview.api.fullscreen()}>
        Test fullscreen
      </button>
      <div style={{ padding: '1rem', display: 'flex', gap: '1rem' }}>
        <button onClick={() => navigate('/features/gvatool')}>
          GVATool
        </button>
        <button onClick={() => navigate('/system/file_manager')}>
          File Manager
        </button>
        <button onClick={() => navigate('/system/macros')}>
          Macros
        </button>
        <button onClick={() => navigate('/resources/manual')}>
          Manual
        </button>
      </div>
    </div>
  )
}
