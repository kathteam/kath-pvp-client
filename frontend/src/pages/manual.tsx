import { JSX } from "react"
import { useNavigate } from "react-router-dom"

export default function Manual(): JSX.Element {
  const navigate = useNavigate()
  
  return (
    <div>
      <h1>Manual Page</h1>
      <p>This is the manual page of our application.</p>
      <button onClick={() => navigate('/index.html')}>
        Back to Dashboard
      </button>    
    </div>
  )
}
