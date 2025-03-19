import { JSX } from "react"
import { useNavigate } from "react-router-dom"

export default function FileManager(): JSX.Element {
  const navigate = useNavigate()
  
  return (
    <div>
      <h1>File Manager Page</h1>
      <p>This is the file manager page of our application.</p>
      <button onClick={() => navigate('/index.html')}>
        Back to Dashboard
      </button>
    </div>
  )
}
